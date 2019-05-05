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

    Light(const Color &ambient, const Color &diffuse, const Color &specular, const Vector3D &v) :
            ambient(ambient), diffuse(diffuse), specular(specular) {
        if (v.is_vector()) vector = Vector3D::normalise(v);
        else vector = v;
    }

    bool isInf() const {
        return vector.is_vector();
    }

    Light &operator*=(const Matrix &matrix) {
        vector *= matrix;
        return *this;
    }

};

struct Lights : public std::vector<Light> {
    Lights &operator*=(const Matrix &matrix) {
        for (auto &light: *this) {
            light *= matrix;
        }
        return *this;
    }
};


#endif //ENGINE_CMAKE_LIGHT_H
