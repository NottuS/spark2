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




void image_size(int& rows, int& cols, Mat& image){
	Mat::SizeType size = image.GetLargestPossibleRegion().GetSize();
	rows = size[1];
	cols = size[0];
}







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



	void draw_Rect(Mat& input, int y, int x, int sy, int sx);
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

void MatSystem::exec_extractContour(Mat& input){
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

			cout << media << endl;

			this->draw_Rect(input, min_y, min_x, dy, dx);
			if ( media > 290 && media < 320 )
				cout << "  5 centavos\n";
			else if ( media > 270 && media < 290 )
				cout << " 10 centavos\n";
			else if ( media > 340 && media < 360 )
				cout << " 25 centavos\n";
			else if ( media > 320 && media < 340 )
				cout << " 50 centavos\n";
			else if ( media > 370 && media < 400 )
				cout << "100 centavos\n";
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
		//cout << rows << " " << cols << endl;


		Mat& borrado = ms.exec_blur(img);

		Mat& border  = ms.exec_threshold(borrado,100,255);


		Mat& menor = ms.exec_erode(border);

		//MatFloat& sobel = ms.exec_sobel(img);

		ms.exec_extractContour(menor);
	
		//Mat& menor = ms.resize(img, 50, 50);


		ms.imwrite(string("tmp/")+getName(argv[i]), menor);
	}
	return 0;
}













