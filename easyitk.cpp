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


#include <stdio.h>


using namespace std;




typedef itk::Image<unsigned char,2> Mat;
typedef itk::Image<float,2>         MatFloat;
typedef itk::CastImageFilter< MatFloat, Mat > CastFilterType;
typedef itk::BinaryBallStructuringElement<unsigned char,2> StructuringElementType;


#define uchar unsigned char
/*
class Mat2 {
	Mat data;


	uchar at(int y, int x){
	}





};*/







class MatSystem {
	itk::ImageFileReader<Mat>::Pointer      reader;
	itk::ImageFileWriter<Mat>::Pointer      writer;
    itk::ImageFileWriter<MatFloat>::Pointer writerfloat;



	itk::SobelEdgeDetectionImageFilter<Mat,MatFloat>::Pointer         sobel;
	itk::CannyEdgeDetectionImageFilter<MatFloat, MatFloat>::Pointer   canny;

	itk::ContourExtractor2DImageFilter<Mat>::Pointer      contour;
	itk::ResampleImageFilter<Mat, Mat>::Pointer           resample;

	itk::BinomialBlurImageFilter<Mat,Mat>::Pointer  blur;
	itk::BinaryThresholdImageFilter <Mat,Mat>::Pointer threshold;

	itk::BinaryErodeImageFilter <Mat,Mat, StructuringElementType>::Pointer erode;


	CastFilterType::Pointer castfilter;


  public:
	MatSystem(){
		reader      = itk::ImageFileReader<Mat>::New();
		writer      = itk::ImageFileWriter<Mat>::New();
		writerfloat = itk::ImageFileWriter<MatFloat>::New();
		sobel       = itk::SobelEdgeDetectionImageFilter <Mat, MatFloat>::New();
		castfilter  = CastFilterType::New();
		//canny       = itk::CannyEdgeDetectionImageFilter<Mat, Mat>::New();
		contour     = itk::ContourExtractor2DImageFilter<Mat>::New();


		resample  = itk::ResampleImageFilter<Mat, Mat>::New();
		blur      = itk::BinomialBlurImageFilter<Mat,Mat>::New();
		threshold = itk::BinaryThresholdImageFilter <Mat,Mat>::New();

		erode     = itk::BinaryErodeImageFilter <Mat,Mat,StructuringElementType>::New();		
	}


	Mat& imread(std::string url);	
	void imwrite(std::string url, Mat& image);
	void imwrite(std::string url, MatFloat& image);

	Mat& resize(Mat& image, int new_rows, int new_cols);


	MatFloat& exec_sobel(Mat& input);
	Mat&      exec_canny(Mat& input);
	void      exec_extractContour(Mat& input);
	Mat&      exec_blur(Mat& input);
	Mat&      exec_threshold(Mat& input, int lower, int upper);
	Mat&      exec_erode(Mat& input);
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


MatFloat& MatSystem::exec_sobel(Mat& input){
	sobel->SetInput(&input);
	sobel->Update();
	return *sobel->GetOutput();
}


Mat& MatSystem::exec_canny(Mat& input){
	//canny->SetInput(&input);
	//canny->Update();
	//return *canny->GetOutput();


	//ms.cannySegmentation->SetInput( reader2->GetOutput() );
	//cannySegmentation->SetFeatureImage( diffusion->GetOutput() );

}


void MatSystem::exec_extractContour(Mat& input){
	contour->SetInput(&input);
	contour->SetContourValue(25);
	contour->Update();

	std::cout << "There are " << contour->GetNumberOfOutputs() << " contours" << std::endl;
	for(unsigned int i = 0; i < contour->GetNumberOfOutputs(); i++){
		//std::cout << "Contour " << i << ": " << std::endl;
		itk::VectorContainer<unsigned int, itk::ContinuousIndex<double, 2> >::Iterator  it = contour->GetOutput(i)->GetVertexList()->Begin();

		while(it != contour->GetOutput(i)->GetVertexList()->End()){
			std::cout << it->Value() << std::endl;
			it++;
		}
		std::cout << std::endl;
	}

	/*Mat::IndexType pt;
    pt[0] = 100;
    pt[1] = 100;
    
	//unsigned char pixelValue = image->GetPixel(pixelIndex);

	input.SetPixel(pt,255);*/








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










void image_size(int& rows, int& cols, Mat& image){
	Mat::SizeType size = image.GetLargestPossibleRegion().GetSize();
	rows = size[1];
	cols = size[0];
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




int main(int argc, char** argv){



	int rows, cols;
	for (int i=1; i<argc; i++){

		cout << argv[i] << endl;


		MatSystem ms;
		Mat& img = ms.imread( argv[i] );


		image_size(rows,cols,img);
		cout << rows << " " << cols << endl;


		Mat& borrado = ms.exec_blur(img);

		Mat& border  = ms.exec_threshold(borrado,100,255);


		Mat& menor = ms.exec_erode(border);

		//MatFloat& sobel = ms.exec_sobel(img);

		//ms.exec_extractContour(menor);
	
		//Mat& menor = ms.resize(img, 50, 50);

























		ms.imwrite(string("tmp/")+getName(argv[i]), menor);
	}
	return 0;
}













