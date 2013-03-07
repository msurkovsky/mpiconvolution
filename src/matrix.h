#ifndef MATRIX_H
#define MATRIX_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

template<typename T>
class Matrix {

    public:
        Matrix() {
            this->width = 0;
            this->height = 0;
            array = new T[0];
        }

        Matrix(unsigned int width, unsigned int height) {
            this->width = width;
            this->height = height;
            array = new T[width * height];
            // fill matrix by zeros;
            memset(array, 0, width * height * sizeof(T));
        }

        Matrix(unsigned int width, unsigned int height,const T *values) {
            this->width = width;
            this->height = height;

            array = new T[width * height];
            memcpy(array, values, width * height * sizeof(T));
        }

        ~Matrix() {
            if (array) {
                delete [] array;
                array = NULL;
            }
        }

        Matrix<T> &operator=(const Matrix<T> &m) {
            if (this != &m) {
                if (width * height != m.width * m.height) {
                    delete [] array;

                    width = m.width;
                    height = m.height;
                    array = new T[width * height];
                }
                memcpy(array, m.array, width * height * sizeof(T));
            }
            return *this;
        }

        T get(unsigned int x, unsigned int y) {
            if (check_coordinates(x, y)) {
                // throw exception
            }

            return array[y * width + x];
        }

        void set(unsigned int x, unsigned int y, T value) {
            if (check_coordinates(x, y)) {
                // throw exception
            }

            array[y * width + x] = value;
        }

        unsigned int get_width() {
            return width;
        }

        unsigned int get_height() {
            return height;
        }

        void clean() {  // fill matrix with zero
            memset(array, 0, width * height * sizeof(T));
        }

        void print(int precision=1) {

            for (int i = 0; i < width; i++) {
                for (int j = 0; j < height; j++) {
                    std::cout.width(precision);
                    std::cout << array[j * width + i] << " ";
                }
                std::cout << std::endl;
            }
        }

        void convolve(Matrix<T> &mask) {

            int mask_width = mask.get_width();
            int mask_height = mask.get_height();

            std::vector<Triple> mask_values;
            for (int i = 0; i < mask_width; i++) {
                for (int j = 0; j < mask_height; j++) {
                    if (mask.get(i, j) != 0) {
                        mask_values.push_back(Triple(i-mask_width/2, j-mask_height/2, mask.get(i, j)));
                    }
                }
            }

            int mask_values_size = mask_values.size();

            int i_start = mask_width / 2;
            int i_end = width - i_start;

            int j_start = mask_height / 2;
            int j_end = height - j_start;

            T *out_data = new T[width * height];
            memset(out_data, 0, width * height * sizeof(T));

            int x, y;
            Triple mask_value;
            for (int i = i_start; i < i_end; i++) {
                for (int j = j_start; j < j_end; j++) {
                    T sum = 0;
                    for (int m = 0; m < mask_values_size; m++) {
                        mask_value = mask_values[m];
                        x = i + mask_value.x;
                        y = j + mask_value.y;
                        sum += (array[y * width + x] * mask_value.value);
                    }
                    out_data[j * width + i] = sum;
                }
            }

            T *tmp = array;
            array = out_data;
            delete [] tmp;
        }

        Matrix<T> *divide(const unsigned int x_divide, const unsigned int y_divide,
                const unsigned int w_overlay, const unsigned int h_overlay) {

            int count = x_divide * y_divide;
            unsigned int original_block_width = width / x_divide;
            unsigned int original_block_height = height / y_divide;
            unsigned int original_block_size = original_block_width * original_block_height;

            unsigned int block_width = original_block_width + 2 * w_overlay;
            unsigned int block_height = original_block_height + 2 * h_overlay;
            unsigned int block_size = block_width * block_height;

            T **blocks;
            blocks = new T*[count];
            for (int t = 0; t < count; t++ ) {
                blocks[t] = new T[block_size];
                memset(blocks[t], 0, sizeof(T) * block_size);
            }

            // fill blocks
            Range x_ranges[x_divide];
            Range y_ranges[y_divide];
            int begin, end;
            for (int i = 0; i < x_divide; i++) {
               begin = i * original_block_width;
               end = begin + original_block_width - 1;
               x_ranges[i] = Range(begin - w_overlay, end + w_overlay);
            }

            for (int i = 0; i < y_divide; i++) {
                begin = i * original_block_height;
                end = begin + original_block_height - 1;
                y_ranges[i] = Range(begin - h_overlay, end + h_overlay);
            }

            for (int y = 0; y < y_divide; y++) {
                Range y_range = y_ranges[y];
                for (int j = y_range.low, j1=0; j <= y_range.heigh; j++, j1++) {

                    if (j < 0 || j >= height) continue;

                    for (int x = 0; x < x_divide; x++) {
                        Range x_range = x_ranges[x];
                        int task = y * x_divide + x;

                        if (x_range.low < 0) {
                            // shift dest pointer
                            int end_overlay = w_overlay;
                            if (x_range.heigh >= width) {
                                // if an overlay happens from the top and the bottom together.
                                end_overlay = 2 * w_overlay;
                            }
                            memcpy(&blocks[task][j1 * block_width + w_overlay],
                                   &array[j * width], // there is + 0 instead of x_range.low
                                   sizeof(T) * (block_width - end_overlay));
                        } else if (x_range.heigh >= width) {
                            // cut block size
                            memcpy(&blocks[task][j1 * block_width],
                                   &array[j * width + x_range.low],
                                   sizeof(T) * (block_width - w_overlay));
                        } else {
                            memcpy(&blocks[task][j1 * block_width],
                                   &array[j * width + x_range.low],
                                   sizeof(T) * block_width);
                        }
                    }
                }
            }

            Matrix<T> *matrices = new Matrix<T>[count];
            for (int i = 0; i < count; i++) {
                Matrix<T>  matrix(block_width, block_height, blocks[i]);
                matrices[i] = matrix;
                delete [] blocks[i];
            }
            delete [] blocks;

            return matrices;
        }

        char * serialize() {
            char *buffer = new char[2*sizeof(unsigned int) + width*height*sizeof(T)];
            memcpy(buffer, &width,  sizeof(unsigned int));
            buffer += sizeof(unsigned int);

            memcpy(buffer, &height, sizeof(unsigned int));
            buffer += sizeof(unsigned int);

            memcpy(buffer, array, width*height*sizeof(T));
        }

        void deserialize(const char *mem) {
            int width, height;
            &width = (unsigned int*) mem;
            mem += sizeof(unsigned int);
            &height = (unsigned int*) mem;
            mem += sizeof(unsigned int);
            std::cout << "Widht: " << width << std::endl;
            std::cout << "Height: " << height << std::endl;
        }

        static Matrix<T> join(const unsigned int width,
                              const unsigned int height,
                              const unsigned int x_divide,
                              const unsigned int y_divide,
                              Matrix<T> * matrices,
                              const unsigned int x_overlay=0,
                              const unsigned int y_overlay=0) {

            int block_width = width / x_divide;
            int block_height = height / y_divide;

            Matrix<T> out(width, height);

            for (int i = 0; i < width; i++) {
                int x = i / block_width;
                for (int j = 0; j < height; j++) {
                    int y = j / block_height;
                    int block_idx = y * x_divide + x;

                    out.set(i, j, matrices[block_idx].get(
                        (i % block_width) + x_overlay,
                        (j % block_height) + y_overlay));
                }
            }

            return out;
        }

    protected:
        struct Triple {
            Triple() { }

            Triple(int x, int y, T value) {
                this->value = value;
                this->x = x;
                this->y = y;
            }

            int x;
            int y;
            T value;
        };

        struct Range {

            int low;
            int heigh;

            Range() {
                this->low = 0;
                this->heigh = 0;
            }

            Range(int low, int heigh) {
                this->low = low;
                this->heigh = heigh;
            }
        };

        unsigned int width;
        unsigned int height;
        T *array;

    private:
        bool check_coordinates(int x, int y) {
            return !(x < 0 || x >= width || y < 0 || y >= height);
        }

        T *get_raw_data() {
            return array;
        }
};

#endif // MATRIX_H
