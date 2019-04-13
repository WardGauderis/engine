#include "easy_image.h"
#include "ini_configuration.h"
#include "Line2D.h"

#include <fstream>

enum render {
	wire, zbuf, triangle
};

img::EasyImage introColorRectangle(const ini::Configuration &configuration) {
	int width = configuration["ImageProperties"]["width"];
	int height = configuration["ImageProperties"]["height"];
	img::EasyImage image(static_cast<unsigned int>(width), static_cast<unsigned int>(height));
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			image(static_cast<unsigned int>(i), static_cast<unsigned int>(j)).red = i;
			image(static_cast<unsigned int>(i), static_cast<unsigned int>(j)).green = j;
			image(static_cast<unsigned int>(i), static_cast<unsigned int>(j)).blue = (i + j) % 256;
		}
	}
	return image;
}

img::EasyImage introBlocks(const ini::Configuration &configuration) {
	int width = configuration["ImageProperties"]["width"];
	int height = configuration["ImageProperties"]["height"];
	img::EasyImage image(static_cast<unsigned int>(width), static_cast<unsigned int>(height));
	int nrXBlocks = configuration["BlockProperties"]["nrXBlocks"];
	int nrYBlocks = configuration["BlockProperties"]["nrYBlocks"];
	std::vector<double> colorWhite = configuration["BlockProperties"]["colorWhite"];
	std::vector<double> colorBlack = configuration["BlockProperties"]["colorBlack"];
	Color white = colorWhite;
	Color black = colorBlack;
	bool invertColors = configuration["BlockProperties"]["invertColors"];
	int blockWidth = static_cast<int>(round(width / nrXBlocks));
	int blockHeight = static_cast<int>(round(height / nrYBlocks));
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; ++y) {
			int xCo = x / blockWidth;
			int yCo = y / blockHeight;
			if ((xCo + yCo + invertColors) % 2 == 0) {
				image(static_cast<unsigned int>(x), static_cast<unsigned int>(y)) = white;
			} else {
				image(static_cast<unsigned int>(x), static_cast<unsigned int>(y)) = black;
			}
		}
	}
	return image;
}

void quadLines(img::EasyImage &image, unsigned int nrLines, unsigned int height, unsigned int width,
			   const Color &line, unsigned int lineY, unsigned int lineX) {
	auto deltaH = static_cast<unsigned int>(round(height / (nrLines - 1.0)));
	auto deltaW = static_cast<unsigned int>(round(width / (nrLines - 1.0)));
	unsigned int x = 0;
	unsigned int y = 0;
	for (unsigned int i = 0; i < nrLines; ++i) {
		image.draw_line(x, lineY, lineX, y, line);
		x += deltaW;
		y += deltaH;
	}
}

img::EasyImage introLines(const ini::Configuration &configuration) {
	int width = configuration["ImageProperties"]["width"];
	int height = configuration["ImageProperties"]["height"];
	img::EasyImage image(static_cast<unsigned int>(width), static_cast<unsigned int>(height));
	std::string figure = configuration["LineProperties"]["figure"];
	int nrLines = configuration["LineProperties"]["nrLines"];
	std::vector<double> backgroundcolor = configuration["LineProperties"]["backgroundcolor"];
	std::vector<double> lineColor = configuration["LineProperties"]["lineColor"];
	Color background = backgroundcolor;
	Color line = lineColor;
	image.clear(background);
	if (figure == "QuarterCircle") {
		quadLines(image, static_cast<unsigned int>(nrLines), static_cast<unsigned int>(height),
				  static_cast<unsigned int>(width), line, static_cast<unsigned int>(height - 1), 0);
	} else if (figure == "Eye") {
		quadLines(image, static_cast<unsigned int>(nrLines), static_cast<unsigned int>(height),
				  static_cast<unsigned int>(width), line, static_cast<unsigned int>(height - 1), 0);
		quadLines(image, static_cast<unsigned int>(nrLines), static_cast<unsigned int>(height),
				  static_cast<unsigned int>(width), line, 0, static_cast<unsigned int>(width - 1));
	} else if (figure == "Diamond") {
		auto deltaH = static_cast<unsigned int>(round(height / (2 * (nrLines - 1.0))));
		auto deltaW = static_cast<unsigned int>(round(width / (2 * (nrLines - 1.0))));
		for (unsigned int i = 0; i < static_cast<unsigned int>(nrLines); ++i) {
			image.draw_line((width / 2) + i * deltaW, static_cast<unsigned int>(height / 2),
							static_cast<unsigned int>(width / 2),
							(height / 2) + ((nrLines - 1 - i) * deltaH), line);
			image.draw_line((width / 2) + i * deltaW, static_cast<unsigned int>(height / 2),
							static_cast<unsigned int>(width / 2),
							(height / 2) - ((nrLines - 1 - i) * deltaH), line);
			image.draw_line((width / 2) - ((nrLines - 1 - i) * deltaW), static_cast<unsigned int>(height / 2),
							static_cast<unsigned int>(width / 2),
							(height / 2) + i * deltaH, line);
			image.draw_line((width / 2) - ((nrLines - 1 - i) * deltaW), static_cast<unsigned int>(height / 2),
							static_cast<unsigned int>(width / 2),
							(height / 2) - i * deltaH, line);
		}
	}
	return image;
}

img::EasyImage lSystem2D(const ini::Configuration &configuration) {
	const std::string inputfile = configuration["2DLSystem"]["inputfile"];
	const int size = configuration["General"]["size"];
	std::vector<double> backgroundcolor = configuration["General"]["backgroundcolor"];
	std::vector<double> linecolor = configuration["2DLSystem"]["color"];
	Color background = backgroundcolor;
	Color line = linecolor;
	LParser::LSystem2D lSystem2D;
	std::ifstream inputStream(inputfile);
	inputStream >> lSystem2D;
	inputStream.close();
	Lines2D lines{lSystem2D, line};
	return lines.draw((unsigned int) size, background, false);
}

img::EasyImage draw3D(const ini::Configuration &configuration, const render type) {
	const int size = configuration["General"]["size"];
	std::vector<double> backgroundcolor = configuration["General"]["backgroundcolor"];
	Color background = backgroundcolor;
	int nrFigures = configuration["General"]["nrFigures"];
	std::vector<double> eyeP = configuration["General"]["eye"];
	Vector3D eye = Vector3D::point(eyeP);
	Figures figures;
	for (int i = nrFigures - 1; i >= 0; --i) {
		std::string figureType = configuration["Figure" + std::to_string(i)]["type"];
		double scale = configuration["Figure" + std::to_string(i)]["scale"];
		double x = M_PI * configuration["Figure" + std::to_string(i)]["rotateX"].as_double_or_default(0) / 180;
		double y = M_PI * configuration["Figure" + std::to_string(i)]["rotateY"].as_double_or_default(0) / 180;
		double z = M_PI * configuration["Figure" + std::to_string(i)]["rotateZ"].as_double_or_default(0) / 180;
		std::vector<double> centerPoint = configuration["Figure" + std::to_string(i)]["center"];
		Vector3D center = Vector3D::point(centerPoint);
		std::vector<double> color = configuration["Figure" + std::to_string(i)]["color"];
		Figure figure;
		if (figureType == "LineDrawing") {
			int nrPoints = configuration["Figure" + std::to_string(i)]["nrPoints"];
			int nrLines = configuration["Figure" + std::to_string(i)]["nrLines"];
			for (int j = 0; j < nrPoints; ++j) {
				std::vector<double> point = configuration["Figure" + std::to_string(i)]["point" + std::to_string(j)];
				figure.addPoint(Vector3D::point(point));
			}
			for (int k = 0; k < nrLines; ++k) {
				std::vector<int> line = configuration["Figure" + std::to_string(i)]["line" + std::to_string(k)];
				figure.addFace(line);
			}
		} else if (figureType == "Cube") {
			figure = Figure::cube();
		} else if (figureType == "Tetrahedron") {
			figure = Figure::tetrahedron();
		} else if (figureType == "Octahedron") {
			figure = Figure::octahedron();
		} else if (figureType == "Icosahedron") {
			figure = Figure::icosahedron();
		} else if (figureType == "Dodecahedron") {
			figure = Figure::dodecahedron();
		} else if (figureType == "Cylinder") {
			const int n = configuration["Figure" + std::to_string(i)]["n"];
			const double height = configuration["Figure" + std::to_string(i)]["height"];
			figure = Figure::cylinder(n, height);
		} else if (figureType == "Cone") {
			const int n = configuration["Figure" + std::to_string(i)]["n"];
			const double height = configuration["Figure" + std::to_string(i)]["height"];
			figure = Figure::cone(n, height);
		} else if (figureType == "Sphere") {
			const int n = configuration["Figure" + std::to_string(i)]["n"];
			figure = Figure::sphere(n);
		} else if (figureType == "Torus") {
			const double R = configuration["Figure" + std::to_string(i)]["R"];
			const double r = configuration["Figure" + std::to_string(i)]["r"];
			const int n = configuration["Figure" + std::to_string(i)]["n"];
			const int m = configuration["Figure" + std::to_string(i)]["m"];
			figure = Figure::Torus(R, r, n, m);
		} else if (figureType == "3DLSystem") {
			const std::string inputfile = configuration["Figure" + std::to_string(i)]["inputfile"];
			LParser::LSystem3D lSystem3D;
			std::ifstream inputStream(inputfile);
			inputStream >> lSystem3D;
			inputStream.close();
			figure = lSystem3D;
		} else if (figureType == "FractalCube") {
			const int nrIterations = configuration["Figure" + std::to_string(i)]["nrIterations"];
			const double fractalScale = configuration["Figure" + std::to_string(i)]["fractalScale"];
			figures += Figures::fractal(Figure::cube(), nrIterations, fractalScale, color);
		} else if (figureType == "FractalTetrahedron") {
			const int nrIterations = configuration["Figure" + std::to_string(i)]["nrIterations"];
			const double fractalScale = configuration["Figure" + std::to_string(i)]["fractalScale"];
			figures += Figures::fractal(Figure::tetrahedron(), nrIterations, fractalScale, color);
		} else if (figureType == "FractalOctahedron") {
			const int nrIterations = configuration["Figure" + std::to_string(i)]["nrIterations"];
			const double fractalScale = configuration["Figure" + std::to_string(i)]["fractalScale"];
			figures += Figures::fractal(Figure::octahedron(), nrIterations, fractalScale, color);
		} else if (figureType == "FractalIcosahedron") {
			const int nrIterations = configuration["Figure" + std::to_string(i)]["nrIterations"];
			const double fractalScale = configuration["Figure" + std::to_string(i)]["fractalScale"];
			figures += Figures::fractal(Figure::icosahedron(), nrIterations, fractalScale, color);
		} else if (figureType == "FractalDodecahedron") {
			const int nrIterations = configuration["Figure" + std::to_string(i)]["nrIterations"];
			const double fractalScale = configuration["Figure" + std::to_string(i)]["fractalScale"];
			figures += Figures::fractal(Figure::dodecahedron(), nrIterations, fractalScale, color);
		}
		if (!figure.getPoints().empty()) {
			figure *= scaleFigure(scale) * rotateX(x) * rotateY(y) * rotateZ(z) * translate(center);
			figure.setColor(color);
			figures.addFigure(figure);
		}
	}
	figures *= eyePoint(eye);
	if (type == wire) {
		Lines2D lines = figures;
		return lines.draw((unsigned int) size, background, false);
	} else if (type == zbuf) {
		Lines2D lines = figures;
		return lines.draw((unsigned int) size, background, true);
	} else if (type == triangle) {
		figures.triangulate();
		return figures.draw((unsigned int) size, background);
	}
	return img::EasyImage();
}

img::EasyImage generate_image(const ini::Configuration &configuration) {
	std::string type = configuration["General"]["type"].as_string_or_die();
	if (type == "IntroColorRectangle") {
		return introColorRectangle(configuration);
	} else if (type == "IntroBlocks") {
		return introBlocks(configuration);
	} else if (type == "IntroLines") {
		return introLines(configuration);
	} else if (type == "2DLSystem") {
		return lSystem2D(configuration);
	} else if (type == "Wireframe") {
		return draw3D(configuration, wire);
	} else if (type == "ZBufferedWireframe") {
		return draw3D(configuration, zbuf);
	} else if (type == "ZBuffering") {
		return draw3D(configuration, triangle);
	}
	return img::EasyImage();
}

int main(int argc, char const *argv[]) {
	int retVal = 0;
	try {
		for (int i = 1; i < argc; ++i) {
			ini::Configuration conf;
			try {
				std::ifstream fin(argv[i]);
				fin >> conf;
				fin.close();
			}
			catch (ini::ParseException &ex) {
				std::cerr << "Error parsing file: " << argv[i] << ": " << ex.what() << std::endl;
				retVal = 1;
				continue;
			}

			img::EasyImage image = generate_image(conf);
			if (image.get_height() > 0 && image.get_width() > 0) {
				std::string fileName(argv[i]);
				std::string::size_type pos = fileName.rfind('.');
				if (pos == std::string::npos) {
					//filename does not contain a '.' --> append a '.bmp' suffix
					fileName += ".bmp";
				} else {
					fileName = fileName.substr(0, pos) + ".bmp";
				}
				try {
					std::ofstream f_out(fileName.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
					f_out << image;

				}
				catch (std::exception &ex) {
					std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
					retVal = 1;
				}
			} else {
				std::cout << "Could not generate image for " << argv[i] << std::endl;
			}
		}
	}
	catch (const std::bad_alloc &exception) {
		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
		//(Unless of course you are already consuming the maximum allowed amount of memory)
		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
		std::cerr << "Error: insufficient memory" << std::endl;
		retVal = 100;
	}
	return retVal;
}
