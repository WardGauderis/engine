#include <utility>

//============================================================================
// @name        : Figure.cpp
// @author      : Ward Gauderis
// @date        : 3/03/2019
//============================================================================

#include <cfloat>
#include "Figure.h"
#include <algorithm>

Matrix scaleFigure(const double scale) {
	Matrix scaler;
	for (int i = 1; i < 4; ++i) {
		scaler(i, i) = scale;
	}
	return scaler;
}

Matrix rotateX(double angle) {
	Matrix rotator;
	for (int i = 2; i < 4; ++i) {
		rotator(i, i) = cos(angle);
	}
	rotator(2, 3) = sin(angle);
	rotator(3, 2) = -sin(angle);
	return rotator;
}

Matrix rotateY(double angle) {
	Matrix rotator;
	for (int i = 1; i < 4; i += 2) {
		rotator(i, i) = cos(angle);
	}
	rotator(1, 3) = -sin(angle);
	rotator(3, 1) = sin(angle);
	return rotator;
}

Matrix rotateZ(double angle) {
	Matrix rotator;
	for (int i = 1; i < 3; i++) {
		rotator(i, i) = cos(angle);
	}
	rotator(1, 2) = sin(angle);
	rotator(2, 1) = -sin(angle);
	return rotator;
}

Matrix translate(const Vector3D &vector) {
	Matrix translator;
	translator(4, 1) = vector.x;
	translator(4, 2) = vector.y;
	translator(4, 3) = vector.z;
	return translator;
}

Matrix eyePoint(const Vector3D &eyepoint) {
	double r = sqrt(pow(eyepoint.x, 2) + pow(eyepoint.y, 2) + pow(eyepoint.z, 2));
	double phi = acos(eyepoint.z / r);
	double theta = atan2(eyepoint.y, eyepoint.x);
	auto trans = Vector3D::point(0, 0, -r);
	Matrix eye = rotateZ((-M_PI / 2) - theta) * rotateX(-phi) * translate(trans);
//	Matrix eye;
//	eye(1, 1) = -sin(theta);
//	eye(1, 2) = -cos(theta) * cos(phi);
//	eye(1, 3) = cos(theta) * sin(phi);
//	eye(2, 1) = cos(theta);
//	eye(2, 2) = -sin(theta) * cos(phi);
//	eye(2, 3) = sin(theta) * sin(phi);
//	eye(3, 2) = sin(phi);
//	eye(3, 3) = cos(phi);
//	eye(4, 3) = -r;
	return eye;
}

Figure &Figure::operator*=(const Matrix &matrix) {
	for (auto &point: points) {
		point *= matrix;
	}
	return *this;
}

const std::vector<Vector3D> &Figure::getPoints() const {
	return points;
}

const std::vector<Face> &Figure::getFaces() const {
	return faces;
}

const Color &Figure::getColor() const {
	return color;
}

void Figure::addPoint(const Vector3D &vector) {
	points.push_back(vector);
}

void Figure::addFace(const Face &face) {
	faces.push_back(face);
}

void Figure::setColor(const Color &newColor) {
	Figure::color = newColor;
}

Figure Figure::cube() {
	Figure figure;
	figure.points.reserve(8);
	figure.addPoint(Vector3D::point(1, -1, -1));
	figure.addPoint(Vector3D::point(-1, 1, -1));
	figure.addPoint(Vector3D::point(1, 1, 1));
	figure.addPoint(Vector3D::point(-1, -1, 1));
	figure.addPoint(Vector3D::point(1, 1, -1));
	figure.addPoint(Vector3D::point(-1, -1, -1));
	figure.addPoint(Vector3D::point(1, -1, 1));
	figure.addPoint(Vector3D::point(-1, 1, 1));
	figure.faces.reserve(6);
	figure.addFace({0, 4, 2, 6});
	figure.addFace({4, 1, 7, 2});
	figure.addFace({1, 5, 3, 7});
	figure.addFace({5, 0, 6, 3});
	figure.addFace({6, 2, 7, 3});
	figure.addFace({0, 5, 1, 4});
	return figure;
}

Figure Figure::tetrahedron() {
	Figure figure;
	figure.addPoint(Vector3D::point(1, -1, -1));
	figure.addPoint(Vector3D::point(-1, 1, -1));
	figure.addPoint(Vector3D::point(1, 1, 1));
	figure.addPoint(Vector3D::point(-1, -1, 1));
	figure.addFace({0, 1, 2});
	figure.addFace({1, 3, 2});
	figure.addFace({0, 3, 1});
	figure.addFace({0, 2, 3});
	return figure;
}

Figure Figure::octahedron() {
	Figure figure;
	figure.points.reserve(6);
	figure.addPoint(Vector3D::point(1, 0, 0));
	figure.addPoint(Vector3D::point(0, 1, 0));
	figure.addPoint(Vector3D::point(-1, 0, 0));
	figure.addPoint(Vector3D::point(0, -1, 0));
	figure.addPoint(Vector3D::point(0, 0, -1));
	figure.addPoint(Vector3D::point(0, 0, 1));
    figure.faces.reserve(8);
	figure.addFace({0, 1, 5});
	figure.addFace({1, 2, 5});
	figure.addFace({2, 3, 5});
	figure.addFace({3, 0, 5});
	figure.addFace({1, 0, 4});
	figure.addFace({2, 1, 4});
	figure.addFace({3, 2, 4});
	figure.addFace({0, 3, 4});
	return figure;
}

Figure Figure::icosahedron() {
	Figure figure;
    figure.points.reserve(12);
	for (int i = 1; i <= 12; ++i) {
		if (i == 1) {
			figure.addPoint(Vector3D::point(0, 0, sqrt(5) / 2));
		} else if (i <= 6) {
			figure.addPoint(Vector3D::point(cos((i - 2) * 2 * M_PI / 5), sin((i - 2) * 2 * M_PI / 5), 0.5));
		} else if (i <= 11) {
			figure.addPoint(
					Vector3D::point(cos(M_PI / 5 + (i - 7) * 2 * M_PI / 5), sin(M_PI / 5 + (i - 7) * 2 * M_PI / 5),
									-0.5));
		} else {
			figure.addPoint(Vector3D::point(0, 0, -sqrt(5) / 2));
		}
	}
    figure.faces.reserve(20);
	figure.addFace({0, 1, 2});
	figure.addFace({0, 2, 3});
	figure.addFace({0, 3, 4});
	figure.addFace({0, 4, 5});
	figure.addFace({0, 5, 1});
	figure.addFace({1, 6, 2});
	figure.addFace({2, 6, 7});
	figure.addFace({2, 7, 3});
	figure.addFace({3, 7, 8});
	figure.addFace({3, 8, 4});
	figure.addFace({4, 8, 9});
	figure.addFace({4, 9, 5});
	figure.addFace({5, 9, 10});
	figure.addFace({5, 10, 1});
	figure.addFace({1, 10, 6});
	figure.addFace({11, 7, 6});
	figure.addFace({11, 8, 7});
	figure.addFace({11, 9, 8});
	figure.addFace({11, 10, 9});
	figure.addFace({11, 6, 10});
	return figure;
}

Figure Figure::buckyball() {
	Figure ico = Figure::icosahedron();
	std::vector<Vector3D> newPoints;
	newPoints.reserve(60);
	std::vector<Face> pentagons;
	for (int i = 0; i < 12; ++i) {
		Face newPent = {};
		int mod = 0;
		std::vector<Face *> adjacent;
		for (auto &face: ico.faces) {
			auto it = find(face.point_indexes.begin(), face.point_indexes.end(), i);
			if (it != face.point_indexes.end()) {
				adjacent.push_back(&face);
			}
		}
		sort(adjacent, i);
		for (auto &face: adjacent) {
			auto it = find(face->point_indexes.begin(), face->point_indexes.end(), i);
			newPent.point_indexes.emplace_back(i + (mod * 12));
			int nextPoint = *((it + 1 == face->point_indexes.end()) ? face->point_indexes.begin() : it + 1);
			newPoints[i + (mod * 12)] =
					(2.0 / 3.0) * ico.getPoints()[i] + (1.0 / 3.0) * ico.getPoints()[nextPoint % 12];
			*it = i + (mod * 12);
			face->point_indexes.insert(it, i + (((++mod) % 5) * 12));
		}
		pentagons.emplace_back(newPent);
	}
	ico.points = newPoints;
//	ico.faces.insert(ico.faces.end(), pentagons.begin(), pentagons.end());
	ico.faces.reserve(pentagons.size());
	std::move(pentagons.begin(), pentagons.end(), std::back_inserter(ico.faces));
	return ico;
}


Figure Figure::dodecahedron() {
	Figure ico = icosahedron();
	Figure figure;
	figure.points.reserve(20);
	for (const auto &face: ico.getFaces()) {
		figure.addPoint((ico.getPoints()[face.point_indexes[0]] + ico.getPoints()[face.point_indexes[1]] +
						 ico.getPoints()[face.point_indexes[2]]) / 3);
	}
    figure.faces.reserve(12);
	figure.addFace({0, 1, 2, 3, 4});
	figure.addFace({0, 5, 6, 7, 1});
	figure.addFace({1, 7, 8, 9, 2});
	figure.addFace({2, 9, 10, 11, 3});
	figure.addFace({3, 11, 12, 13, 4});
	figure.addFace({4, 13, 14, 5, 0});
	figure.addFace({19, 18, 17, 16, 15});
	figure.addFace({19, 14, 13, 12, 18});
	figure.addFace({18, 12, 11, 10, 17});
	figure.addFace({17, 10, 9, 8, 16});
	figure.addFace({16, 8, 7, 6, 15});
	figure.addFace({15, 6, 5, 14, 19});
	return figure;
}

Figure Figure::cylinder(const int n, const double height) {
	Figure figure;
    figure.points.reserve(2*n);
	for (int i = 0; i < n; ++i) {
		figure.addPoint(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), 0));
	}
	for (int i = 0; i < n; ++i) {
		figure.addPoint(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), height));
	}
    figure.faces.reserve(n + 2);
	for (int j = 0; j < n; ++j) {
		figure.addFace({j, (j + 1) % n, ((j + 1) % n) + n, j + n});
	}
	std::vector<int> points1;
    points1.reserve(n);
	for (int k = 0; k < n; ++k) {
		points1.push_back(n - 1 - k);
	}
	figure.addFace(points1);
	std::vector<int> points2;
	points2.reserve(n);
	for (int k = 0; k < n; ++k) {
		points2.push_back(2 * n - 1 - k);
	}
	figure.addFace(points2);
	return figure;
}

Figure Figure::cone(const int n, const double height) {
	Figure figure;
    figure.points.reserve(n + 1);
	for (int i = 0; i < n; ++i) {
		figure.addPoint(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), 0));
	}
	figure.addPoint(Vector3D::point(0, 0, height));
    figure.faces.reserve(n);
	for (int j = 0; j < n; ++j) {
		figure.addFace({j, (j + 1) % n, n});
	}
	std::vector<int> points;
    points.reserve(n);
	points.reserve((unsigned long) n);
	for (int k = 0; k < n; ++k) {
		points.push_back(n - 1 - k);
	}
	figure.addFace(points);
	return figure;
}

Figure Figure::sphere(const int n) {
	Figure figure = icosahedron();
	for (int j = 0; j < n; ++j) {
        figure.points.reserve(2*figure.points.size());
        figure.faces.reserve(3*figure.faces.size());
		const int size = static_cast<const int>(figure.getFaces().size());
		for (int i = 0; i < size; i++) {
			const Face face = figure.getFaces()[i];
			const Vector3D A = figure.getPoints()[face.point_indexes[0]];
			const Vector3D B = figure.getPoints()[face.point_indexes[1]];
			const Vector3D C = figure.getPoints()[face.point_indexes[2]];
			const Vector3D D = (A + B) / 2;
			const Vector3D E = (A + C) / 2;
			const Vector3D F = (B + C) / 2;
			figure.addPoint(D);
			figure.addPoint(E);
			figure.addPoint(F);
			figure.faces[i] = {face.point_indexes[0], static_cast<int>(figure.getPoints().size() - 3),
                               static_cast<int>(figure.getPoints().size() - 2)};
			figure.addFace({face.point_indexes[1], static_cast<int>(figure.getPoints().size() - 1),
							static_cast<int>(figure.getPoints().size() - 3)});
			figure.addFace({face.point_indexes[2], static_cast<int>(figure.getPoints().size() - 2),
							static_cast<int>(figure.getPoints().size() - 1)});
			figure.addFace(
					{static_cast<int>(figure.getPoints().size() - 3), static_cast<int>(figure.getPoints().size() - 1),
					 static_cast<int>(figure.getPoints().size() - 2)});
		}
	}
	figure.normalize();
	return figure;
}

void Figure::normalize() {
	for (auto &point: points) {
		point.normalise();
	}
}

Figure Figure::Torus(const double R, const double r, const int n, const int m) {
	Figure figure;
    figure.points.reserve(n * m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			const double u = 2 * i * M_PI / n;
			const double v = 2 * j * M_PI / m;
			const double x = (R + r * cos(v)) * cos(u);
			const double y = (R + r * cos(v)) * sin(u);
			const double z = r * sin(v);
			figure.addPoint(Vector3D::point(x, y, z));
		}
	}
    figure.faces.reserve(n*m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			figure.addFace(
					{i * m + j, ((i + 1) % n) * m + j, ((i + 1) % n) * m + ((j + 1) % m), i * m + ((j + 1) % m)});
		}
	}
	return figure;
}

Figure::Figure() = default;

Figure::Figure(const LParser::LSystem3D &lSystem) {
	Vector3D start = Vector3D::point(0, 0, 0);
	addPoint(start);
	int prevPoint = 0;
	Vector3D H = Vector3D::vector(1, 0, 0);
	Vector3D L = Vector3D::vector(0, 1, 0);
	Vector3D U = Vector3D::vector(0, 0, 1);
	unsigned int nr = lSystem.get_nr_iterations();
	const double angle = lSystem.get_angle() * M_PI / 180;
	std::stack<std::tuple<Vector3D, Vector3D, Vector3D, Vector3D, int>> brackets;
	for (const auto &symbol:lSystem.get_initiator()) {
		drawCharacter(nr, start, H, L, U, lSystem, symbol,
					  brackets, angle, prevPoint);
	}
}

void Figure::drawCharacter(unsigned int nr, Vector3D &start, Vector3D &H, Vector3D &L, Vector3D &U,
						   const LParser::LSystem3D &lSystem, const char character,
						   std::stack<std::tuple<Vector3D, Vector3D, Vector3D, Vector3D, int>> &brackets,
						   const double angle, int &prevPoint) {
	Vector3D old;
	switch (character) {
		case '+':
			old = H;
			H = H * cos(angle) + L * sin(angle);
			L = -old * sin(angle) + L * cos(angle);
			break;
		case '-':
			old = H;
			H = H * cos(-angle) + L * sin(-angle);
			L = -old * sin(-angle) + L * cos(-angle);
			break;
		case '^':
			old = H;
			H = H * cos(angle) + U * sin(angle);
			U = -old * sin(angle) + U * cos(angle);
			break;
		case '&':
			old = H;
			H = H * cos(-angle) + U * sin(-angle);
			U = -old * sin(-angle) + U * cos(-angle);
			break;
		case '\\':
			old = L;
			L = L * cos(angle) - U * sin(angle);
			U = old * sin(angle) + U * cos(angle);
			break;
		case '/':
			old = L;
			L = L * cos(-angle) - U * sin(-angle);
			U = old * sin(-angle) + U * cos(-angle);
			break;
		case '(':
			brackets.emplace(start, H, L, U, prevPoint);
			break;
		case ')':
			start = std::get<0>(brackets.top());
			H = std::get<1>(brackets.top());
			L = std::get<2>(brackets.top());
			U = std::get<3>(brackets.top());
			prevPoint = std::get<4>(brackets.top());
			brackets.pop();
			break;
		case '|':
			H = -H;
			L = -L;
			break;
		default:
			if (nr == 0) {
				start += H;
				addPoint(start);
				if (lSystem.draw(character)) {
					addFace({prevPoint, static_cast<int>(getPoints().size() - 1)});
				}
				prevPoint = static_cast<int>(getPoints().size() - 1);
			} else {
				for (const auto &symbol:lSystem.get_replacement(character)) {
					drawCharacter(nr - 1, start, H, L, U, lSystem, symbol, brackets, angle, prevPoint);
				}
			}
			break;
	}
}

void Figure::triangulate() {
	std::vector<Face> backup = std::move(faces);
	for (auto &face: backup) {
		face.triangulate(faces);
	}
}

Figure Figure::operator*(const Matrix &matrix) const {
	Figure fig = *this;
	fig *= matrix;
	return fig;
}

void Figure::sort(std::vector<Face *> &faces, int index) {
	for (int face = 0; face < faces.size() - 2; ++face) {
		int point = 0;
		while (point + 1 < faces[face]->point_indexes.size() && faces[face]->point_indexes[point + 1] != index) point++;
		int prev = faces[face]->point_indexes[point];
		int toSwap = 0;
		for (int test = face + 1; test < faces.size(); ++test) {
			if (std::find(faces[test]->point_indexes.begin(), faces[test]->point_indexes.end(), prev) !=
				faces[test]->point_indexes.end()) {
				toSwap = test;
				break;
			}
		}
		std::swap(faces[face + 1], faces[toSwap]);
	}
}

////Figure Figure::mengerSponge(const int iter) {
////	Figures menger;
////	menger.addFigure(Figure::cube() * rotateX(1) * rotateY(1));
////	double scale = 1;
////	for (int i = 0; i < iter; ++i) {
////		Figures newMenger;
////		for (const auto &cube: menger.getFigures()) {
////			Figure newCube = cube * scaleFigure(1.0/3);
////			std::vector<Vector3D> points = {
////					cube.getPoints()[5] + Vector3D::vector(0, 0, 0),
////					cube.getPoints()[5] + Vector3D::vector(scale, 0, 0),
////					cube.getPoints()[5] + Vector3D::vector(2 * scale, 0, 0),
////					cube.getPoints()[5] + Vector3D::vector(0, scale, 0),
////					cube.getPoints()[5] + Vector3D::vector(2 * scale, scale, 0),
////					cube.getPoints()[5] + Vector3D::vector(0, 2 * scale, 0),
////					cube.getPoints()[5] + Vector3D::vector(scale, 2 * scale, 0),
////					cube.getPoints()[5] + Vector3D::vector(2 * scale, 2 * scale, 0),
////
////					cube.getPoints()[5] + Vector3D::vector(0, 0, scale),
////					cube.getPoints()[5] + Vector3D::vector(2 * scale, 0, scale),
////					cube.getPoints()[5] + Vector3D::vector(0, 2 * scale, scale),
////					cube.getPoints()[5] + Vector3D::vector(2 * scale, 2 * scale, scale),
////
////					cube.getPoints()[5] + Vector3D::vector(0, 0, 2 * scale),
////					cube.getPoints()[5] + Vector3D::vector(scale, 0, 2 * scale),
////					cube.getPoints()[5] + Vector3D::vector(2 * scale, 0, 2 * scale),
////					cube.getPoints()[5] + Vector3D::vector(0, scale, 2*scale),
////					cube.getPoints()[5] + Vector3D::vector(2 * scale, scale, 2*scale),
////					cube.getPoints()[5] + Vector3D::vector(0, 2 * scale, 2*scale),
////					cube.getPoints()[5] + Vector3D::vector(scale, 2 * scale, 2*scale),
////					cube.getPoints()[5] + Vector3D::vector(2 * scale, 2 * scale, 2*scale),
////			};
////			for (const auto &point: points) {
////				newMenger.addFigure(newCube * translate(cube.getPoints()[5] - point));
////			}
////		}
////		menger = std::move(newMenger);
////		scale /= 3;
////	}
////	Figure mengerSimple = Figure(menger, false);
////	int i = mengerSimple.faces.size();
////	return mengerSimple;
////}
//
//Figure::Figure(Figures &figures, bool removeDoubleFaces) {
//	for (auto &figure: figures.figures) {
//		faces.reserve(figure.faces.size());
//		for (auto &face: figure.faces) {
//			for (auto &index: face.point_indexes) {
//				index += points.size();
//			}
//			faces.push_back(std::move(face));
//		}
//		points.reserve(figure.points.size());
//		std::move(figure.points.begin(), figure.points.end(), std::back_inserter(points));
//	}
//	if (removeDoubleFaces) {
//		Vector3D toFind;
//		auto compare = [&toFind](const Vector3D &a) {
//			return (a.x == toFind.x) && (a.y == toFind.y) && (a.z == toFind.z);
//		};
//		for (int i = 0; i < points.size(); i++) {
//			if (points[i].x != NAN) {
//				toFind = points[i];
//				auto it = std::find_if(points.begin() + i + 1, points.end(), compare);
//				while (it != points.end()) {
//					*it = Vector3D::point({NAN, NAN, NAN});
//					int pos = it - points.begin();
//					for (auto &face: faces) {
//						auto it2 = find(face.point_indexes.begin(), face.point_indexes.end(), pos);
//						if (it2 != face.point_indexes.end()) {
//							*it2 = i;
//						}
//					}
//					it = std::find_if(it + 1, points.end(), compare);
//				}
//			}
//		}
//		for (auto it = faces.begin(), next = it; it != faces.end(); it = next) {
//			next++;
//			bool erased = false;
//			for (auto it2 = it + 1, next2 = it2; it2 != faces.end(); it2 = next2) {
//				next2++;
//				int count = 0;
//				for (auto &point: (*it).point_indexes) {
//					if (std::find((*it2).point_indexes.begin(), (*it2).point_indexes.end(), point) !=
//						(*it2).point_indexes.end()) {
//						count++;
//					}
//				}
//				if (count > 2) {
//					next2 = faces.erase(it2);
//					erased = true;
//				}
//			}
//			if (erased) {
//				next = faces.erase(it);
//			}
//		}
//	}
//	color = figures.getFigures().front().getColor();
//}

img::EasyImage Figures::draw(unsigned int size, const Color &background) const {
	auto xMax = -DBL_MAX;
	auto xMin = DBL_MAX;
	auto yMax = -DBL_MAX;
	auto yMin = DBL_MAX;

	for (const auto &figure: figures) {
		for (const auto &point: figure.getPoints()) {
			const double x = (point.x) / (-point.z);
			const double y = (point.y) / (-point.z);
			if (x > xMax) xMax = x;
			if (x < xMin) xMin = x;
			if (y > yMax) yMax = y;
			if (y < yMin) yMin = y;
		}
	}
	const double xRange = xMax - xMin;
	const double yRange = yMax - yMin;
	const double xImage = size * (xRange / (std::max(xRange, yRange)));
	const double yImage = size * (yRange / (std::max(xRange, yRange)));
	const double d = 0.95 * xImage / xRange;
	const double DCx = d * ((xMin + xMax) / 2);
	const double DCy = d * ((yMin + yMax) / 2);
	const double dx = (xImage / 2) - DCx;
	const double dy = (yImage / 2) - DCy;
	img::EasyImage image(static_cast<unsigned int>(round(xImage)), static_cast<unsigned int>(round(yImage)));
	ZBuffer buffer(image.get_width(), image.get_height());
	image.clear(background);
	for (const auto &figure: figures) {
		for (const auto &triangle: figure.getFaces()) {
			image.draw_triangle(buffer, figure.getPoints()[triangle.point_indexes[0]],
								figure.getPoints()[triangle.point_indexes[1]],
								figure.getPoints()[triangle.point_indexes[2]], d, dx, dy, figure.getColor());
		}
	}
	return image;
}

Figures &Figures::operator*=(const Matrix &matrix) {
	for (auto &figure: figures) {
		figure *= matrix;
	}
	return *this;
}

const std::forward_list<Figure> &Figures::getFigures() const {
	return figures;
}

void Figures::addFigure(Figure &&figure) {
	figures.emplace_front(std::move(figure));
}

void Figures::triangulate() {
	for (auto &figure: figures) {
		figure.triangulate();
	}
}

Figures Figures::fractal(Figure &figure, const int iter, const double scale) {
	Figures figs;
	figs.addFigure(std::move(figure));
	for (int i = 0; i < iter; i++) {
		std::forward_list<Figure> newIt;
		for (auto &fig : figs.figures) {
			Figure temp = fig * scaleFigure(1 / scale);
			for (int j = 0; j < fig.getPoints().size(); j++) {
				newIt.emplace_front(temp * translate(fig.getPoints()[j] - temp.getPoints()[j]));
			}
		}
		figs.figures = std::move(newIt);
	}
	return figs;
}

Figures &Figures::operator+=(Figures &&figs) {
	figures.splice_after(figures.before_begin(), figs.figures);
	return *this;
}

Figures Figures::operator*(const Matrix &matrix) const {
	Figures fig = *this;
	fig *= matrix;
	return fig;
}

Figures Figures::mengerSponge(const int iter) {
	Figures mengerSponge;
	mengerRec(mengerSponge, iter, 1, Vector3D::point(-1, -1, -1));
	return mengerSponge;
}

void Figures::mengerRec(Figures &figs, const int iter, double scale, const Vector3D &corner) {
	if (iter == 0) {
		figs.figures.push_front(Figure::cube() * scaleFigure(scale) * translate(corner-Vector3D::point(-1, -1, -1)));
	} else {
		scale /= 3;
		double s = scale * 2;
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(0, 0, 0));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(s, 0, 0));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(2 * s, 0, 0));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(0, s, 0));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(2 * s, s, 0));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(0, 2 * s, 0));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(s, 2 * s, 0));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(2 * s, 2 * s, 0));

		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(0, 0, s));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(2 * s, 0, s));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(0, 2 * s, s));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(2 * s, 2 * s, s));

		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(0, 0, 2 * s));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(s, 0, 2 * s));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(2 * s, 0, 2 * s));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(0, s, 2 * s));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(2 * s, s, 2 * s));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(0, 2 * s, 2 * s));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(s, 2 * s, 2 * s));
		mengerRec(figs, iter - 1, scale, corner + Vector3D::vector(2 * s, 2 * s, 2 * s));
	}
}

void Figures::setColor(const Color &newColor) {
	for (auto &figure: figures) {
		figure.setColor(newColor);
	}
}

Face::Face(
		const std::initializer_list<int> &point_indexes) : point_indexes(point_indexes) {}

Face::Face(std::vector<int>
		   point_indexes) : point_indexes(std::move(point_indexes)) {}

void Face::triangulate(std::vector<Face> &faces) {
	faces.reserve(point_indexes.size() - 2);
	for (unsigned int i = 1; i < point_indexes.size() - 1; ++i) {
		faces.emplace_back(std::initializer_list<int>{point_indexes[0], point_indexes[i], point_indexes[i + 1]});
	}
}
