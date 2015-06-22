#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSobelEdgeDetectionImageFilter.h>
#include <itkCannyEdgeDetectionImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkContourExtractor2DImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkBinomialBlurImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkLineIterator.h>
#include <itkLineConstIterator.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkNeighborhoodOperatorImageFilter.h>
#include <itkSobelOperator.h>
#include <itkAdaptiveHistogramEqualizationImageFilter.h>


#include <itkOtsuThresholdImageFilter.h>



#include <stdio.h>
#include <iostream>
#include <vector>


using namespace std;




typedef itk::Image<unsigned char,2> Mat;
typedef itk::Image<float,2>         MatFloat;
typedef itk::CastImageFilter< MatFloat, Mat > CastFilterType;
typedef itk::BinaryBallStructuringElement<unsigned char,2> StructuringElementType;
typedef itk::NeighborhoodOperatorImageFilter<Mat, MatFloat> NeighborhoodOperatorImageFilterType;
	

#define uchar unsigned char

/*class Mat2 {
	Mat data;


	void set(int y, int x, char value){
		Mat::IndexType pt;
		pt[0] = x;
		pt[1] = y;
		input.SetPixel(pt,value);
	}
};*/


void image_set(Mat& input, int y, int x, uchar value){
	Mat::IndexType pt;
	pt[0] = x;
	pt[1] = y;
	input.SetPixel(pt,value);
}


float image_get(MatFloat& input, int y, int x){
	Mat::IndexType pt;
	pt[0] = x;
	pt[1] = y;
	return input.GetPixel(pt);
}



void image_size(int& rows, int& cols, Mat& image){
	Mat::SizeType size = image.GetLargestPossibleRegion().GetSize();
	rows = size[1];
	cols = size[0];
}



class Rect{
  public:
	int x,y,sx,sy;

	Rect(int y,int x,int sy,int sx){
		this->x = x;
		this->y = y;
		this->sx = sx;
		this->sy = sy;
	}
};





class MatSystem {
	itk::ImageFileReader<Mat>::Pointer      reader;
	itk::ImageFileWriter<Mat>::Pointer      writer;
    itk::ImageFileWriter<MatFloat>::Pointer writerfloat;



	itk::SobelEdgeDetectionImageFilter<Mat,MatFloat>::Pointer         sobel;

	NeighborhoodOperatorImageFilterType::Pointer sobel1;
	NeighborhoodOperatorImageFilterType::Pointer sobel2;

	itk::CannyEdgeDetectionImageFilter<MatFloat, MatFloat>::Pointer   canny;

	itk::ContourExtractor2DImageFilter<Mat>::Pointer      contour;
	itk::ResampleImageFilter<Mat, Mat>::Pointer           resample;

	itk::BinomialBlurImageFilter<Mat,Mat>::Pointer  blur;
	itk::BinaryThresholdImageFilter <Mat,Mat>::Pointer threshold;

	itk::OtsuThresholdImageFilter<Mat,Mat>::Pointer otsu;

	itk::BinaryErodeImageFilter <Mat,Mat, StructuringElementType>::Pointer erode;


	itk::RegionOfInterestImageFilter<Mat,Mat>::Pointer roi;

	itk::AdaptiveHistogramEqualizationImageFilter<Mat>::Pointer equalize;

	CastFilterType::Pointer castfilter;


  public:
	MatSystem(){
		reader      = itk::ImageFileReader<Mat>::New();
		writer      = itk::ImageFileWriter<Mat>::New();
		writerfloat = itk::ImageFileWriter<MatFloat>::New();
		sobel       = itk::SobelEdgeDetectionImageFilter <Mat, MatFloat>::New();
		castfilter  = CastFilterType::New();
		canny       = itk::CannyEdgeDetectionImageFilter<MatFloat, MatFloat>::New();
		contour     = itk::ContourExtractor2DImageFilter<Mat>::New();


		resample  = itk::ResampleImageFilter<Mat, Mat>::New();
		blur      = itk::BinomialBlurImageFilter<Mat,Mat>::New();
		threshold = itk::BinaryThresholdImageFilter <Mat,Mat>::New();

		erode     = itk::BinaryErodeImageFilter <Mat,Mat,StructuringElementType>::New();

		roi       = itk::RegionOfInterestImageFilter<Mat,Mat>::New();

		sobel1    = NeighborhoodOperatorImageFilterType::New();
		sobel2    = NeighborhoodOperatorImageFilterType::New();

		otsu      = itk::OtsuThresholdImageFilter <Mat,Mat>::New();

		equalize  = itk::AdaptiveHistogramEqualizationImageFilter<Mat>::New();
	}


	Mat& imread(std::string url);	
	void imwrite(std::string url, Mat& image);
	void imwrite(std::string url, MatFloat& image);

	Mat& resize(Mat& image, int new_rows, int new_cols);


	MatFloat& exec_sobel(Mat& input);
	MatFloat& exec_canny(MatFloat& input);
	void      exec_extractContour(vector<Rect>& output, Mat& input);
	Mat&      exec_blur(Mat& input);
	Mat&      exec_threshold(Mat& input, int lower, int upper);
	Mat&      exec_erode(Mat& input);

	Mat&      exec_thresholdOtsu(Mat& input);

	Mat&      get_roi(Mat& input, Rect rect);

	void draw_Rect(Mat& input, int y, int x, int sy, int sx);

	Mat& exec_equalize(Mat& input);

	MatFloat& exec_sobelX(Mat& input);
	MatFloat& exec_sobelY(Mat& input);
};



Mat&  MatSystem::imread(std::string url){
	reader->SetFileName( url );
	reader->Update();
	return *reader->GetOutput();
}

void MatSystem::imwrite(std::string url, Mat& image){
	writer->SetFileName(url);
	writer->SetInput(  &image  );
	writer->Update();
}


void MatSystem::imwrite(std::string url, MatFloat& image){
	castfilter->SetInput(&image);

	writer->SetFileName(url);
	writer->SetInput( castfilter->GetOutput() );
	writer->Update();
}


Mat& MatSystem::resize(Mat& image, int new_rows, int new_cols){
	Mat::SizeType outputSize;
	outputSize[0] = new_cols;
	outputSize[1] = new_rows;

	Mat::SpacingType outputSpacing;
	Mat::SizeType inputSize = image.GetLargestPossibleRegion().GetSize();
	outputSpacing[0] = image.GetSpacing()[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
	outputSpacing[1] = image.GetSpacing()[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));

	resample->SetInput(&image);
	resample->SetSize(outputSize);
	resample->SetOutputSpacing(outputSpacing);
	resample->SetTransform(  itk::IdentityTransform<double, 2>::New()  );
	resample->UpdateLargestPossibleRegion();

	return *resample->GetOutput();
}


MatFloat& MatSystem::exec_sobel(Mat& input ){
	sobel->SetInput(&input);
	sobel->Update();
	return *sobel->GetOutput();
}


MatFloat& MatSystem::exec_canny(MatFloat& input){
	canny->SetInput(&input);
	canny->Update();
	return *canny->GetOutput();


	//ms.cannySegmentation->SetInput( reader2->GetOutput() );
	//cannySegmentation->SetFeatureImage( diffusion->GetOutput() );

}


void MatSystem::draw_Rect(Mat& input, int y, int x, int sy, int sx){
	for (int i=0; i<sx; i++){
		image_set(input,y,x+i,64);
	}
	for (int i=0; i<sx; i++){
		image_set(input,y+sy,x+i,64);
	}
	for (int i=0; i<sy; i++){
		image_set(input,y+i,x,64);
	}
	for (int i=0; i<sy; i++){
		image_set(input,y+i,x+sx,64);
	}
}

void MatSystem::exec_extractContour(vector<Rect>& output, Mat& input){
	output.clear();
	contour->SetInput(&input);
	contour->SetContourValue(25);
	contour->Update();

	//std::cout << "There are " << contour->GetNumberOfOutputs() << " contours" << std::endl;
	for(unsigned int i = 0; i < contour->GetNumberOfOutputs(); i++){
		//std::cout << "Contour " << i << ": " << std::endl;
		itk::VectorContainer<unsigned int, itk::ContinuousIndex<double, 2> >::Iterator  it = contour->GetOutput(i)->GetVertexList()->Begin();

		if ( contour->GetOutput(i)->GetVertexList()->size() > 0 ){
			itk::ContinuousIndex<double, 2u> ini = it->Value();
			int min_x = ini[0];
			int max_x = ini[0];
			int min_y = ini[1];
			int max_y = ini[1];
			while(it != contour->GetOutput(i)->GetVertexList()->End()){
				itk::ContinuousIndex<double, 2u> pt = it->Value();
				if ( pt[0] < min_x )
					min_x = pt[0];
				if ( pt[0] > max_x )
					max_x = pt[0];
				if ( pt[1] < min_y )
					min_y = pt[1];
				if ( pt[1] > max_y )
					max_y = pt[1];
				it++;
			}

			int dx = max_x - min_x;
			int dy = max_y - min_y;
			int media = ( dx + dy )/2;


			//cout << dy << " " << dx << " " << dy*dy+dx*dx << " ";
			if ( abs(dy-dx) > 20 || media < 200 ){
				continue;
				//cout << "quadrado";
			}
			//cout << endl;

			//cout << media << endl;

			output.push_back( Rect(min_y,min_x,dy,dx) );
			//this->draw_Rect(input, min_y, min_x, dy, dx);
			if ( media >= 300 && media <= 325 )
				cout << "  5 centavos - ";
			else if ( media > 263 && media <= 299 )
				cout << " 10 centavos - ";
			else if ( media > 340 && media <=376 )
				cout << " 25 centavos - ";
			else if ( media>=326 && media <= 340 )
				cout << " 50 centavos - ";
			else if ( media >= 377 && media < 410 )
				cout << "100 centavos - ";
			cout << media << endl;
		}

	}



	//return contour->GetOutput();
}


Mat& MatSystem::exec_blur(Mat& input){
	blur->SetInput(&input);
	blur->Update();
	return *blur->GetOutput();
}


Mat& MatSystem::exec_threshold(Mat& input, int lower, int upper){
	threshold->SetInput(&input);
	threshold->SetLowerThreshold(lower);
	threshold->SetUpperThreshold(upper);
	threshold->SetInsideValue(255);
	threshold->SetOutsideValue(0);
	threshold->Update();

	return *threshold->GetOutput();
}



Mat& MatSystem::exec_erode(Mat& input){
	StructuringElementType structuringelement;
	structuringelement.SetRadius(2);
	structuringelement.CreateStructuringElement();

	erode->SetInput(&input);
	erode->SetKernel(structuringelement);
	erode->Update();
	return *erode->GetOutput();
}


Mat&  MatSystem::exec_thresholdOtsu(Mat& input){
	otsu->SetInput(&input);
	otsu->Update();
	return *otsu->GetOutput();
}



Mat& MatSystem::get_roi(Mat& input, Rect rect){
	Mat::IndexType start;
	start[0] = rect.x;
	start[1] = rect.y;

	Mat::SizeType size;
	size[0] = rect.sx;
	size[1] = rect.sy;

	Mat::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	roi->SetRegionOfInterest(desiredRegion);
	roi->SetInput(&input);
	roi->Update();

	return *roi->GetOutput();
}




typedef itk::SobelOperator<float, 2> SobelOperatorType;

MatFloat& MatSystem::exec_sobelX(Mat& input){
	SobelOperatorType sobelOperator;
	itk::Size<2> radius;
	radius.Fill(1);
	sobelOperator.SetDirection(0);
	sobelOperator.CreateToRadius(radius);
	sobel1->SetOperator(sobelOperator);
	sobel1->SetInput(&input);
	sobel1->Update();
	return *sobel1->GetOutput();
}


MatFloat& MatSystem::exec_sobelY(Mat& input){
	SobelOperatorType sobelOperator;
	itk::Size<2> radius;
	radius.Fill(1);
	sobelOperator.SetDirection(1);
	sobelOperator.CreateToRadius(radius);
	sobel2->SetOperator(sobelOperator);
	sobel2->SetInput(&input);
	sobel2->Update();
	return *sobel2->GetOutput();
}







string getName(string url){
	string buffer;
	for (int i=0; i<url.size(); i++){
		char c = url[i];
		if ( c == '/' ){
			buffer = "";
		} else
			buffer += c;
	}
	return buffer;
}


Mat& MatSystem::exec_equalize(Mat& input){
	equalize->SetInput(&input);
	equalize->SetRadius(1);
	equalize->Update();
	return *equalize->GetOutput();
}



int main(int argc, char** argv){

	system("mkdir -p tmp");

	int rows, cols;
	for (int i=1; i<argc; i++){

		cout << argv[i] << endl;


		MatSystem ms;
		Mat& img = ms.imread( argv[i] );


		image_size(rows,cols,img);
		//cout << rows << " " << cols << endl;


		Mat& borrado = ms.exec_blur(img);

		//Mat& canelson = ms.exec_canny(borrado);


		//Mat& equalized = ms.exec_equalize(borrado);
		Mat& border  = ms.exec_threshold(borrado,120,255);

		//Mat& border  = ms.exec_thresholdOtsu(borrado);

		
		//Mat& menor = ms.exec_erode(border);





		vector<Rect> rectpkg;
		ms.exec_extractContour(rectpkg, border);
		for (int j=0; j<rectpkg.size(); j++){
			ms.draw_Rect(border, rectpkg[j].y, rectpkg[j].x, rectpkg[j].sy, rectpkg[j].sx);
		}
		

		ms.imwrite(string("tmp/")+getName(argv[i]), border);
		//Mat& menor = ms.resize(img, 50, 50);


		//ms.imwrite(string("tmp/")+getName(argv[i]), menor);
	}
	return 0;
}






/*for (int j=0; j<rectpkg.size(); j++){
			Mat& coin = ms.get_roi(img,rectpkg[j]);
			ms.imwrite(getName("a.jpg"), coin);

			MatFloat& coin_border_x = ms.exec_sobelX(coin);
			MatFloat& coin_border_y = ms.exec_sobelY(coin);

			int border_rows,border_cols;
			image_size(border_rows,border_cols,coin);


			int counter[8];
			for (int i=0; i<8; i++)
				counter[i] = 0;

			double t_tx=0.0, t_ty=0.0;

			for (int iy=0; iy<border_rows; iy++){
				for (int ix=0; ix<border_cols; ix++){
					float th = 0;
					float ay = image_get(coin_border_y,iy,ix);
					float ax = image_get(coin_border_x,iy,ix);
					th = atan2(ay,ax);
					int direction =   ( (th*180.0/3.141592)/45 ) + 3;
					counter[direction] += fabs(ay) + fabs(ax);





					if ( fabs(ax) + fabs(ay) > 100 ){
						t_tx += ax;
						t_ty += ay;
					}
					
				}
			}

			double desc[8];
			double total = border_cols*border_rows;
			for (int i=0; i<8; i++){
				desc[i] = counter[i]/total;
				cout << desc[i] << endl;
			}

		}*/







