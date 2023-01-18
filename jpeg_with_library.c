#include "library/stb/stb_image.h"
#include "library/stb/stb_image_write.h"

int main() {
    int width, height, channels;
    unsigned char *image = stbi_load("naruto.png", &width, &height, &channels, 0);

    if (!image) {
        printf("Failed to load image\n");
        return 1;
    }

    if (!stbi_write_jpg("naruto.jpg", width, height, channels, image, 100)) {
        printf("Failed to write JPEG\n");
        return 1;
    }

    stbi_image_free(image);
    return 0;
}
