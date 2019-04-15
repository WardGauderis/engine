//============================================================================
// @name        : Figure.h
// @author      : Ward Gauderis
// @date        : 3/03/2019
//============================================================================

#ifndef ENGINE_CMAKE_FIGURE_H
#define ENGINE_CMAKE_FIGURE_H

#include <vector>
#include <forward_list>
#include "l_parser/l_parser.h"
#include "vector/vector3d.h"
#include "easy_image.h"
#include "ini_configuration.h"
#include <stack>

Matrix scaleFigure(double scale);

Matrix rotateX(double angle);

Matrix rotateY(double angle);

Matrix rotateZ(double angle);

Matrix translate(const Vector3D &vector);

Matrix eyePoint(const Vector3D &eyepoint);

struct Face {
	std::vector<int> point_indexes;

	Face(const std::initializer_list<int> &point_indexes);

	Face(std::vector<int> point_indexes);

	void triangulate(std::vector<Face> &faces);
};

class Figures;

class Figure {
	std::vector<Vector3D> points;
	std::vector<Face> faces;
	Color color;

	static void sort(std::vector<Face *> &faces, int index);

public:
	Figure &operator*=(const Matrix &matrix);

	Figure operator*(const Matrix &matrix) const;

	const std::vector<Vector3D> &getPoints() const;

	const std::vector<Face> &getFaces() const;

	void setColor(const Color &newColor);

	const Color &getColor() const;

	void addPoint(const Vector3D &vector);

	void addFace(const Face &face);

	void deleteFace(int i);

	void normalize();

	static Figure cube();

	static Figure tetrahedron();

	static Figure octahedron();

	static Figure icosahedron();

	static Figure buckyball();

	static Figure dodecahedron();

	static Figure cylinder(int n, double height);

	static Figure cone(int n, double height);

	static Figure sphere(int n);

	static Figure Torus(double R, double r, int n, int m);

//	static Figure mengerSponge(int iter);

	Figure(const LParser::LSystem3D &lSystem);

	void triangulate();

	Figure();

//	explicit Figure(Figures& figures, bool removeDoubleFaces);

	void drawCharacter(unsigned int nr, Vector3D &start, Vector3D &H, Vector3D &L, Vector3D &U,
					   const LParser::LSystem3D &lSystem, char character,
					   std::stack<std::tuple<Vector3D, Vector3D, Vector3D, Vector3D, int>> &brackets,
					   double angle,
					   int &prevPoint);
};

class Figures {
	std::forward_list<Figure> figures;

	static void mengerRec(Figures &figs, const int iter, const double scale, const Vector3D &corner);
public:
	Figures &operator*=(const Matrix &matrix);

	Figures &operator+=(Figures &&figs);

	Figures operator*(const Matrix &matrix) const;

	const std::forward_list<Figure> &getFigures() const;

	void triangulate();

	void addFigure(Figure &&figure);

	img::EasyImage draw(unsigned int size, const Color &background) const;

	static Figures fractal(Figure &figure, int iter, double scale);

	static Figures mengerSponge(int iter);

	void setColor(const Color &newColor);

//	friend Figure::Figure(Figures& figures, bool removeDoubleFaces);
};

#endif //ENGINE_CMAKE_FIGURE_H
