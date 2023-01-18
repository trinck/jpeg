#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265
#define MAX_SYMBOLS 256
#define BLOCK_SIZE 8


// Function to read PNG image and convert it to a 2D array of pixels
unsigned char** readPNG(char* fileName, int* width, int* height, int numColors, int resolution) {
    // Open the PNG file
    FILE* pngFile = fopen(fileName, "rb");
    if (!pngFile) {
        printf("Error: Unable to open file %s\n", fileName);
        return NULL;
    }

    // Read the width and height of the image from the PNG file
    fseek(pngFile, 16, SEEK_SET);
    fread(width, 4, 1, pngFile);
    fread(height, 4, 1, pngFile);

    // Read the number of bits per sample
    unsigned char bitsPerSample;
    fseek(pngFile, 25, SEEK_SET);
    fread(&bitsPerSample, 1, 1, pngFile);

    // Calculate the number of colors
    numColors = (int)pow(2, bitsPerSample);


    // Read the resolution
    fseek(pngFile, 38, SEEK_SET);
    fread(&resolution, 4, 1, pngFile);


    // Allocate memory for the image data
    unsigned char** pixels = (unsigned char**)malloc(*height * sizeof(unsigned char*));
    for (int i = 0; i < height; i++) {
        pixels[i] = (unsigned char*)malloc(*width * sizeof(unsigned char));
    }

    // Read the image data into the 2D array
    fseek(pngFile, 41, SEEK_SET);
    for (int i = 0; i < height; i++) {
        fread(pixels[i], width, 1, pngFile);
    }

    // Close the PNG file
    fclose(pngFile);

    return pixels;
}


unsigned char** divideIntoBlocks(unsigned char** pixels, int width, int height, unsigned char*** blocks) {
    // Copy the pixels into the blocks
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int blockX = x / BLOCK_SIZE;
            int blockY = y / BLOCK_SIZE;
            int pixelX = x % BLOCK_SIZE;
            int pixelY = y % BLOCK_SIZE;
            blocks[blockY][blockX][pixelY * BLOCK_SIZE + pixelX] = pixels[y][x];
        }
    }
    return blocks;
}


// Function to perform DCT on a 8x8 block of pixels
double** dct(unsigned char** block, double** dctCoeffs) {
    // Perform DCT on the block of pixels
    for (int u = 0; u < 8; u++) {
        for (int v = 0; v < 8; v++) {
            double cu = (u == 0) ? 1 / sqrt(2) : 1;
            double cv = (v == 0) ? 1 / sqrt(2) : 1;
            double sum = 0;
            for (int x = 0; x < 8; x++) {
                for (int y = 0; y < 8; y++) {
                    double val = (double)block[x][y];
                    double x_val = (double)(cos((2 * x + 1) * u * PI / 16));
                    double y_val = (double)(cos((2 * y + 1) * v * PI / 16));
                    sum += val * x_val * y_val;
                }
            }
            dctCoeffs[u][v] = (cu * cv / 4) * sum;
        }
    }
    return dctCoeffs;
}


// Function to perform quantization on the DCT coefficients
int** quantize(double** dctCoeffs, int** quantCoeffs, int quality) {
    // Define the quantization table
    int quantTable[8][8] = {{16, 11, 10, 16, 24, 40, 51, 61},
                           {12, 12, 14, 19, 26, 58, 60, 55},
                           {14, 13, 16, 24, 40, 57, 69, 56},
                           {14, 17, 22, 29, 51, 87, 80, 62},
                           {18, 22, 37, 56, 68, 109, 103, 77},
                           {24, 35, 55, 64, 81, 104, 113, 92},
                           {49, 64, 78, 87, 103, 121, 120, 101},
                           {72, 92, 95, 98, 112, 100, 103, 99}};

    // Scale the quantization table based on the quality factor
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            quantTable[i][j] = (quantTable[i][j] * quality + 50) / 100;
            if (quantTable[i][j] == 0) {
                quantTable[i][j] = 1;
            }
        }
    }

    // Perform quantization on the DCT coefficients
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            quantCoeffs[i][j] = (int)(dctCoeffs[i][j] / quantTable[i][j]);
        }
    }
    return quantCoeffs;
}



int* zigzag(int** quantCoeffs, int* zigzagCoeffs) {
    // Define the Zigzag pattern
    int pattern[64] = {0, 1, 5, 6, 14, 15, 27, 28,
                       2, 4, 7, 13, 16, 26, 29, 42,
                       3, 8, 12, 17, 25, 30, 41, 43,
                       9, 11, 18, 24, 31, 40, 44, 53,
                       10, 19, 23, 32, 39, 45, 52, 54,
                       20, 22, 33, 38, 46, 51, 55, 60,
                       21, 34, 37, 47, 50, 56, 59, 61,
                       35, 36, 48, 49, 57, 58, 62, 63};

    // Follow the Zigzag pattern and store the DCT coefficients in the 1D array
    for (int i = 0; i < 64; i++) {
        int x = pattern[i] / 8;
        int y = pattern[i] % 8;
        zigzagCoeffs[i] = quantCoeffs[x][y];
    }
    return zigzagCoeffs;
}


void runLengthEncode(double* zigzagCoeffs, int size, int* encodedCoeffs, int encodedSize) {
    int count = 0;
    int zeros = 0;
    for (int i = 0; i < size; i++) {
        if (zigzagCoeffs[i] == 0) {
            zeros++;
        } else {
            if (zeros > 0) {
                encodedCoeffs[count] = 0;
                encodedCoeffs[count + 1] = zeros;
                count += 2;
                zeros = 0;
            }
            encodedCoeffs[count] = zigzagCoeffs[i];
            count++;
        }
    }
    if (zeros > 0) {
        encodedCoeffs[count] = 0;
        encodedCoeffs[count + 1] = zeros;
        count += 2;
    }
    encodedSize = count;
}


typedef struct {
    int symbol;
    int count;
} Symbol;

typedef struct _HuffNode {
    int symbol;
    int count;
    struct _HuffNode* left;
    struct _HuffNode* right;
} HuffNode;

typedef struct {
    HuffNode* root;
    int* codes[MAX_SYMBOLS];
    int codeLengths[MAX_SYMBOLS];
} HuffmanTree;

struct PriorityQueue {
    int size;
    int capacity;
    HuffNode** array;
};


struct bitStream {
    FILE* file;
    unsigned char buffer;
    int index;
};


void initBitStream(struct bitStream* stream, FILE* file) {
    stream->file = file;
    stream->buffer = 0;
    stream->index = 0;
}



void huffmanEncode(int* encodedCoeffs, int encodedSize, struct bitStream* stream) {
    // Create a symbol array to store the counts of each value in the coeffs array
    Symbol symbols[MAX_SYMBOLS];
    for (int i = 0; i < MAX_SYMBOLS; i++) {
        symbols[i].symbol = i;
        symbols[i].count = 0;
    }
    // Count the occurrences of each value in the coeffs array
    for (int i = 0; i < encodedSize; i++) {
        symbols[encodedCoeffs[i]].count++;
    }
    // Create the Huffman tree using the symbol counts
    HuffNode* root = createHuffmanTree(symbols);
    // Generate the Huffman codes for each symbol
    int code[256];
    int codes[MAX_SYMBOLS][256];
    int codeLengths[MAX_SYMBOLS];
    generateCodes(root, code, 0, codes, codeLengths);
    // Write the Huffman codes to the output file using the bit stream
    for (int i = 0; i < encodedSize; i++) {
        int symbol = encodedCoeffs[i];
        for (int j = 0; j < codeLengths[symbol]; j++) {
            writeBit(stream, codes[symbol][j]);
        }
    }
    // Clean up memory used for the Huffman tree
    freeHuffmanTree(root);
}


// Function to write the bit
void writeBit(struct bitStream* stream, int bit) {
    // Add the bit to the buffer
    stream->buffer = (stream->buffer << 1) | bit;
    stream->index++;
    // Flush the buffer to the file if it is full
    if (stream->index == 8) {
        fputc(stream->buffer, stream->file);
        stream->buffer = 0;
        stream->index = 0;
    }
}

HuffNode* createHuffmanTree(Symbol* symbols) {
    // Sort the symbols by count
    qsort(symbols, MAX_SYMBOLS, sizeof(Symbol), compareSymbols);
    // Create a priority queue to hold the Huffman tree nodes
    struct PriorityQueue* queue = createPriorityQueue();

    // Create the leaf nodes and add them to the queue
    for (int i = 0; i < MAX_SYMBOLS; i++) {
        if (symbols[i].count > 0) {
            HuffNode* node = (HuffNode*)malloc(sizeof(HuffNode));
            node->symbol = symbols[i].symbol;
            node->count = symbols[i].count;
            node->left = NULL;
            node->right = NULL;
            enqueue(queue, node);
        }
    }

    // Build the Huffman tree by repeatedly dequeuing the two smallest nodes and
    // creating a new parent node with their counts as the sum
    while (queue->size > 1) {
        HuffNode* left = dequeue(queue);
        HuffNode* right = dequeue(queue);
        HuffNode* parent = (HuffNode*)malloc(sizeof(HuffNode));
        parent->count = left->count + right->count;
        parent->left = left;
        parent->right = right;
        enqueue(queue, parent);
    }
    // The remaining node in the queue is the root of the Huffman tree
    HuffNode* root = dequeue(queue);

    // Clean up
    destroyPriorityQueue(queue);

    return root;
}

// Function to compare the counts of two symbols for sorting
int compareSymbols(const void* a, const void* b) {
    return ((Symbol*)a)->count - ((Symbol*)b)->count;
}

// Function to create a priority queue
struct PriorityQueue* createPriorityQueue() {
    struct PriorityQueue* queue = (struct PriorityQueue*)malloc(sizeof(struct PriorityQueue));
    queue->size = 0;
    queue->capacity = MAX_SYMBOLS;
    queue->array = (HuffNode**)malloc(queue->capacity * sizeof(HuffNode*));
    return queue;
}


// Function to enqueue a node into the priority queue
void enqueue(struct PriorityQueue* queue, HuffNode* node) {
    if (queue->size == queue->capacity) {
        printf("Error: Priority queue is full\n");
        return;
    }
    queue->array[queue->size] = node;
    int current = queue->size;
    int parent = (current - 1) / 2;
    while (parent >= 0 && queue->array[current]->count < queue->array[parent]->count) {
        HuffNode* temp = queue->array[current];
        queue->array[current] = queue->array[parent];
        queue->array[parent] = temp;
        current = parent;
        parent = (current - 1) / 2;
    }
    queue->size++;
}


// Function to dequeue the node with the smallest count
HuffNode* dequeue(struct PriorityQueue* queue) {
    if (queue->size == 0) {
        printf("Error: Priority queue is empty\n");
        return NULL;
    }
    HuffNode* min = queue->array[0];
    queue->array[0] = queue->array[queue->size - 1];
    queue->size--;
    heapify(queue, 0);
    return min;
}


// Function to maintain the heap property of the priority queue
void heapify(struct PriorityQueue* queue, int index) {
    int left = 2 * index + 1;
    int right = 2 * index + 2;
    int smallest = index;
    if (left < queue->size && queue->array[left]->count < queue->array[smallest]->count) {
        smallest = left;    
    }
    if (right < queue->size && queue->array[right]->count < queue->array[smallest]->count) {
        smallest = right;
    }
    if (smallest != index) {
        HuffNode* temp = queue->array[index];
        queue->array[index] = queue->array[smallest];
        queue->array[smallest] = temp;
        heapify(queue, smallest);
    }
}


// Function to destroy the priority queue
void destroyPriorityQueue(struct PriorityQueue* queue) {
    free(queue->array);
    free(queue);
}

void generateCodes(HuffNode* node, int* code, int codeLength, int codes[MAX_SYMBOLS][256], int codeLengths[MAX_SYMBOLS]) {
    if (node->left == NULL && node->right == NULL) {
        // Leaf node, store the code
        memcpy(codes[node->symbol], code, codeLength * sizeof(int));
        codeLengths[node->symbol] = codeLength;
    } else {
        // Non-leaf node, recursively generate codes for the left and right subtrees
        if (node->left != NULL) {
            code[codeLength] = 0;
            generateCodes(node->left, code, codeLength + 1, codes, codeLengths);
        }
        if (node->right != NULL) {
            code[codeLength] = 1;
            generateCodes(node->right, code, codeLength + 1, codes, codeLengths);
        }
    }
}


// Free the Huffman node
void freeHuffmanTree(HuffNode* node) {
    if (node != NULL) {
        freeHuffmanTree(node->left);
        freeHuffmanTree(node->right);
        free(node);
    }
}

// Function to write the JPGHeader 
void writeJPGHeader(FILE* jpgFile, int* encodedCoeffs, int encodedSize, int* zigzagCoeffs, int width, int height, int quality, int resolution, int numColors) {

    // Define the quantization table
    int quantTable[64] = {16, 11, 10, 16, 24, 40, 51, 61,
                            12, 12, 14, 19, 26, 58, 60, 55,
                            14, 13, 16, 24, 40, 57, 69, 56,
                            14, 17, 22, 29, 51, 87, 80, 62,
                            18, 22, 37, 56, 68, 109, 103, 77,
                            24, 35, 55, 64, 81, 104, 113, 92,
                            49, 64, 78, 87, 103, 121, 120, 101,
                            72, 92, 95, 98, 112, 100, 103, 99};


    // Write the SOI marker
    fputc(0xFF, jpgFile);
    fputc(0xD8, jpgFile);

    // Write the JFIF marker
    fputc(0xFF, jpgFile);
    fputc(0xE0, jpgFile);
    fputc(0x00, jpgFile);
    fputc(0x10, jpgFile);
    fwrite("JFIF", 5, 1, jpgFile);

    // Write the resolution
    fputc(0x01, jpgFile);
    fputc(0x02, jpgFile);
    fputc((resolution >> 8) & 0xFF, jpgFile);
    fputc(resolution & 0xFF, jpgFile);
    fputc(0x01, jpgFile);
    fputc(0x02, jpgFile);

    // Write the image size
    fputc(0xFF,jpgFile);
    fputc(0xC0, jpgFile);
    fputc(0x00, jpgFile);
    fputc(0x11, jpgFile);
    fputc(0x08, jpgFile);
    fputc((height >> 8) & 0xFF, jpgFile);
    fputc(height & 0xFF, jpgFile);
    fputc((width >> 8) & 0xFF, jpgFile);
    fputc(width & 0xFF, jpgFile);
    fputc(0x03, jpgFile);
    fputc(numColors & 0xFF, jpgFile);
    fputc(0x00, jpgFile);
    fputc(0x00, jpgFile);
    fputc(0x00, jpgFile);
    fputc(0x00, jpgFile);

    // Write the DQT marker
    fputc(0xFF, jpgFile);
    fputc(0xDB, jpgFile);
    fputc(0x00, jpgFile);
    fputc(0x43, jpgFile);
    fputc(0x00, jpgFile);

    // Write the quantization table
    for (int i = 0; i < 64; i++) {
        int quantValue = (int)(quantTable[i] * quality / 100.0);
        fputc(quantValue & 0xFF, jpgFile);
    }

    // Write the DHT marker
    fputc(0xFF, jpgFile);
    fputc(0xC4, jpgFile);


    // Write the SOS marker
    fputc(0xFF, jpgFile);
    fputc(0xDA, jpgFile);

    // Write the length of the SOS marker
    fputc(0x00, jpgFile);
    fputc(0x0C, jpgFile);

    // Write the number of image components
    fputc(0x03, jpgFile);

    // Write the information about each image component
    /*for (int i = 0; i < 3; i++) { 
        fputc(componentId[i], jpgFile); 
        fputc(componentDcHuffTable[i], jpgFile); 
        fputc(componentAcHuffTable[i], jpgFile); 
    }*/

    // Write length of the Huffman table
    fputc(0x00, jpgFile); 
    fputc((MAX_SYMBOLS >> 8) & 0xFF, jpgFile);
    fputc(MAX_SYMBOLS & 0xFF, jpgFile);

    // Write the information about the Huffman Table
    fputc(0x00, jpgFile);
    fputc(0xAC, jpgFile); 


    Symbol symbols[MAX_SYMBOLS];
    for (int i = 0; i < MAX_SYMBOLS; i++) {
        symbols[i].symbol = i;
        symbols[i].count = 0;
    }
    // Count the occurrences of each value in the coeffs array
    for (int i = 0; i < encodedSize; i++) {
        symbols[encodedCoeffs[i]].count++;
    }
    // Create the Huffman tree using the symbol counts
    HuffNode* root = createHuffmanTree(symbols);
    // Generate the Huffman codes for each symbol
    int code[256];
    int codes[MAX_SYMBOLS][256];
    int codeLengths[MAX_SYMBOLS];
    generateCodes(root, code, 0, codes, codeLengths);

    // Write the number of codes for each bit length
    for (int i = 0; i < encodedSize; i++) { 
        fputc((codeLengths[encodedCoeffs[i]] >> 8) & 0xFF, jpgFile);
        fputc(codeLengths[encodedCoeffs[i]] & 0xFF, jpgFile);
    }

    // Write the Huffman codes for each bit length
    for (int i = 0; i < encodedSize; i++) {
        for (int j = 0; j < codeLengths[encodedCoeffs[i]]; j++) {
            fputc((codes[encodedCoeffs[i]][j] >> 8) & 0xFF, jpgFile);
            fputc(codes[encodedCoeffs[i]][j] & 0xFF, jpgFile);
        }
    }

    freeHuffmanTree(root);

    // Write the end-of-scan marker
    fputc(0x00, jpgFile);
    fputc(0x3F, jpgFile);
    fputc(0x00, jpgFile);

    // Write the spectral selection and approximation data
    fputc(0x00, jpgFile); 
    fputc(0x3F, jpgFile); 
    fputc(0x00, jpgFile);


    // Write the actual scan data
    for (int i = 0; i < 64; i++) { 
        fputc(zigzagCoeffs[i], jpgFile); 
    }

    // Write the length of the scan data
    fputc(0x00, jpgFile);
    fputc((64 >> 8) & 0xFF, jpgFile);
    fputc(64 & 0xFF, jpgFile);


    // Write the EOF marker
    fputc(0x00, jpgFile);
    fputc(0x0C, jpgFile);


}


// Function to write the JPG image using the quantized DCT coefficients
void writeJPG(char* fileName, int* zigzagCoeffs, int* encodedCoeffs, int width, int height, int quality, int resolution, int numColors, int encodedSize) {
    // Open the JPG file
    FILE* jpgFile = fopen(fileName, "wb");

    if (!jpgFile) {
        printf("Error: Unable to open file %s\n", fileName);
        return;
    }
    
    writeJPGHeader(jpgFile, encodedCoeffs, &encodedSize, zigzagCoeffs, &width, &height, &quality, &resolution, &numColors);

    printf("JPG Header");

    // Create a new bit stream to store the encoded data
    struct bitStream* stream;
    initBitStream(stream, jpgFile);

    // Perform Huffman encoding on the encodedCoeffs coefficients
    huffmanEncode(encodedCoeffs, &encodedSize, stream);

    printf("Huffman encode");

    // Write the Huffman encoded data to the JPG file
    if(stream->index>0){
        fwrite(stream->buffer, sizeof(unsigned char), 1, stream->file);
    }
    // Close the JPG file
    fclose(jpgFile);
}





int main(int argc, char* argv[]) {

    if (argc < 3) {
        printf("Error: Please provide a file name\n");
        return 1;
    }
    
    // Read the PNG image and convert it to a 2D array of pixels
    char* inputFile = argv[1];
    int width, height, numColors, resolution;



    unsigned char** pixels = readPNG(inputFile, &width, &height, &numColors, &resolution);
    if (!pixels) {
        return 1;
    }
    for (int i = 0; i < height; i++) {
        pixels[i] = (unsigned char*)malloc(width * sizeof(unsigned char));
    }

    printf("PIXELS");

    // Perform 8*8 blocks
    // Calculate the number of blocks in the image
    int numBlocksX = (int)ceil((double)width / BLOCK_SIZE);
    int numBlocksY = (int)ceil((double)height / BLOCK_SIZE);

    // Allocate memory for the blocks
    unsigned char*** blocks = (unsigned char***)malloc(numBlocksY * sizeof(unsigned char**));
    for (int i = 0; i < numBlocksY; i++) {
        blocks[i] = (unsigned char**)malloc(numBlocksX * sizeof(unsigned char*));
        for (int j = 0; j < numBlocksX; j++) {
            blocks[i][j] = (unsigned char*)malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(unsigned char));
        }
    }
    blocks = divideIntoBlocks(pixels, &width, &height, blocks);

    printf("BLOCS");

    // Perform DCT on the image
    double** dctCoeffs = (double**)malloc(8 * sizeof(double*));
    for (int i = 0; i < 8; i++) {
        dctCoeffs[i] = (double*)malloc(8 * sizeof(double));
    }
    dctCoeffs = dct(blocks, dctCoeffs);

    printf("DCT");

    // Perform quantization on the DCT coefficients
    int quality = 75;
    int** quantCoeffs = (int**)malloc(8 * sizeof(int*));
    for (int i = 0; i < 8; i++) {
        quantCoeffs[i] = (int*)malloc(8 * sizeof(int));
    }
    quantCoeffs = quantize(dctCoeffs, quantCoeffs, &quality);

    printf("QUANTIZATION");

    // Perform the Zigzag pattern on the quantized DCT coefficients
    int* zigzagCoeffs = (int*) malloc(64 * sizeof(int));
    zigzagCoeffs = zigzag(quantCoeffs, zigzagCoeffs);

    printf("ZIGZAG");

    // Perform the runLengthEncode pattern on the quantized DCT coefficients
    int encodedSize;
    int* encodedCoeffs = (int*) malloc(64 * sizeof(int));
    runLengthEncode(zigzagCoeffs, 64, encodedCoeffs, &encodedSize);

    printf("RLE");

    // Write the JPG image
    char* outputFile = argv[2];
    writeJPG(outputFile, zigzagCoeffs, encodedCoeffs, &width, &height, &quality, &resolution, &numColors, &encodedSize);


    free(pixels);

    free(dctCoeffs);

    free(quantCoeffs);

    free(zigzagCoeffs);

    free(inputFile);

    free(outputFile);

    return 0;
}

