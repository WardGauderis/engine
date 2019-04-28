//============================================================================
// @name        : Light.h
// @author      : Ward Gauderis
// @date        : 4/28/19
// @version     : 
// @copyright   : BA1 Informatica - Ward Gauderis - University of Antwerp
// @description : 
//============================================================================
#ifndef ENGINE_CMAKE_LIGHT_H
#define ENGINE_CMAKE_LIGHT_H

#include <forward_list>
#include "Color.h"

struct Light {
    Color ambient;
    Color diffuse;
    Color specular;
    Vector3D vector;

    Light(const Color &ambient, const Color &diffuse, const Color &specular, Vector3D v) :
            ambient(ambient), diffuse(diffuse), specular(specular) {
        if (v.is_vector()) v.normalise();
        vector = v;
    }

    bool isInf() const {
        return vector.is_vector();
    }

    Light &operator*=(const Matrix &matrix) {
        vector *= matrix;
    }

};

struct Lights : public std::forward_list<Light> {
    Lights &operator*=(const Matrix &matrix) {
        for (auto &light: *this) {
            light *= matrix;
        }
    }
};


#endif //ENGINE_CMAKE_LIGHT_H
