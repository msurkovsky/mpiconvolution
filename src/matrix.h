#ifndef MATRIX_H
#define MATRIX_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

template<typename T>
class Matrix {
    
    public:
        Matrix(size_t width, size_t height) {

            this->width = width;
            this->height = height;
            array = new T[width * height];
            // fill matrix by zeros;
            memset(array, 0, width * height * sizeof(T));
        }

        Matrix(size_t width, size_t height, const T *values) {
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

        T get(size_t x, size_t y) {
            if (check_coordinates(x, y)) {
                // throw exception
            }

            return array[y * width + x];
        }

        void set(size_t x, size_t y, T value) {
            if (check_coordinates(x, y)) {
                // throw exception
            }

            array[y * width + x] = value;
        }

        size_t get_width() { return width; }
        size_t get_height() {return height; }

        void clean() {  // fill matrix with zero
            memset(array, 0, width * height * sizeof(T));
        }

        void print() {
            for (int i = 0; i < width; i++) {
                for (int j = 0; j < height; j++) {
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

            int j_start = mask_width / 2;
            int j_end = width - j_start;

            T *out_data = new T[width * height];

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
        size_t width;
        size_t height;
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
