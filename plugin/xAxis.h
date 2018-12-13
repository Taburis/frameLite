

#ifndef xAxis_H
#define xAxis_H

class xAxis{
		public : xAxis (){
				 };
				 ~xAxis(){
						 delete bins;
				 }
				 xAxis (int nbin, float *bin);
				 //	 	 xAxis (int nbin, const float *bin);
				 xAxis (int nbin, float xmin, float xmax);
				 void setAxis(int n, float *bin){
						 bins = new float[n+1];
						 nbin = n;
						 for(int i=0; i<n+1; ++i) bins[i]=bin[i];
				 }
				 int findBin(float x);
				 template <typename T> long BinarySearch(long n, const T *array, T value)
				 {
						 //  Binary search in an array of n values to locate value.
						 //
						 // Array is supposed  to be sorted prior to this call.
						 // If match is found, function returns position of element.
						 // If no match found, function gives nearest element smaller than value.

						 const T* pind;
						 pind = std::lower_bound(array, array + n, value);
						 if ( (pind != array + n) && (*pind == value) )
								 return (pind - array);
						 else
								 return (pind - array - 1);
				 }

		public : 
				 float * bins;
				 int     nbin;
				 float   xmax;
				 float   xmin;
				 char*   name;
};

xAxis::xAxis(int n, float *bin){
		bins = new float[n+1];
		nbin = n;
		for(int i=0; i<n+1; ++i) bins[i]=bin[i];
}
//xAxis::xAxis(int n, const float *bin){
//	bins = new float[n+1];
//	nbin = n;
//	for(int i=0; i<n+1; ++i) bins[i]=bin[i];
//}

xAxis::xAxis(int n, float xmin, float xmax){
		float width =(xmax-xmin)/n;
		bins = new float[n+1];
		nbin = n;
		for(int i=0; i<n+1; ++i) bins[i]=xmin+i*width;
}

int xAxis::findBin(float x){
		return BinarySearch(nbin+1, bins, float(x));
}

#endif

