#ifndef MAT_H
#define MAT_H

#include <QList>
#include <iostream>
using namespace std;

template <int M, int N, typename T = float>
class Mat
{
public:
    Mat();
    Mat(T* data);
    ~Mat();
    T* data();
    T at(int i, int j);
    void set(int i, int j, T v);
    T determinant();
    Mat<N,M,T> transposed();
    Mat<M,N,T> inverse();
    Mat<M-1,N-1,T> cofactor(int i, int j);
    void printData();
    void printData(T* data, int m, int n);
    void printData(T* data, int dim);

private:
    T* m_data;
    T* minorData(T* olddata, int m, int n, int dim);
    T determinantHelper(T* data, int dim);
};

template <int M, int N, typename T>
Mat<M,N,T> operator + (Mat<M,N,T> &m1, Mat<M,N,T> &m2)
{
    T* data = (T*)malloc(M*N*sizeof(T));
    for(int i = 0; i < M*N; i++)
        data[i] = (m1.data())[i]+(m2.data())[i];
    return (Mat<M,N,T>(data));
}

template <int M, int N, typename T>
Mat<M,N,T> operator - (Mat<M,N,T> &m1, Mat<M,N,T> &m2)
{
    T* data = (T*)malloc(M*N*sizeof(T));
    for(int i = 0; i < M*N; i++)
        data[i] = (m1.data())[i]-(m2.data())[i];
    return (Mat<M,N,T>(data));
}

template <int M, int N, typename T>
Mat<M,N,T> operator - (Mat<M,N,T> &m)
{
    T* data = (T*)malloc(M*N*sizeof(T));
    for(int i = 0; i < M*N; i++)
        data[i] = -(m.data())[i];
    return Mat<M,N,T>(data);
}

template <int M, int K, int N, typename T>
Mat<M,N,T> operator * (Mat<M,K,T> m1, Mat<K,N,T> m2)
{
    T* data = (T*)malloc(M*N*sizeof(T));
    for(int i = 0; i < M; i++)
        for(int j = 0; j < N; j++)
        {
            T v = 0;
            T* d1 = m1.data();
            T* d2 = m2.data();
            for(int k = 0; k < K; k++)
            {
                v += d1[i*K+k]*d2[k*N+j];
            }
            data[i*N+j] = v;
        }
    return (Mat<M,N,T>(data));
}

template <int M, int N, typename T>
Mat<M,N,T> operator * (Mat<M,N,T> &m, T v)
{
    T* data = (T*)malloc(M*N*sizeof(T));
    for(int i = 0; i < M*N; i++)
        data[i] = (m.data())[i]*v;
    return (Mat<M,N,T>(data));
}

template <int M, int N, typename T>
Mat<M,N,T> operator * (T v, Mat<M,N,T> &m)
{
    T* data = (T*)malloc(M*N*sizeof(T));
    for(int i = 0; i < M*N; i++)
        data[i] = (m.data())[i]*v;
    return (Mat<M,N,T>(data));
}

template <int M, int N, typename T>
Mat<M,N,T> operator / (Mat<M,N,T> &m, T v)
{
    T* data = (T*)malloc(M*N*sizeof(T));
    for(int i = 0; i < M*N; i++)
        data[i] = (m.data())[i]/v;
    return (Mat<M,N,T>(data));
}

template <int M, int N, typename T>
void Mat<M,N,T>::printData()
{
    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < N; j++)
            cout<<m_data[i*N+j]<<'\t';
        cout<<endl;
    }
    cout<<endl;
}

template <int M, int N, typename T>
void Mat<M,N,T>::printData(T *data, int m, int n)
{
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
            cout<<data[i*n+j]<<'\t';
        cout<<endl;
    }
    cout<<endl;
}

template <int M, int N, typename T>
void Mat<M,N,T>::printData(T *data, int dim)
{
    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
            cout<<data[i*dim+j]<<'\t';
        cout<<endl;
    }
    cout<<endl;
}

template <int M, int N, typename T>
Mat<M,N,T>::Mat()
{
    if(M<1 || N<1)
        throw "Error: M or N can't be less than 1.";
    //m_data = (T*)calloc(M*N,sizeof(T));
    m_data = new T[M*N];
}

template <int M, int N, typename T>
Mat<M,N,T>::Mat(T *data)
{
    if(M<1 || N<1)
        throw "Error: M or N can't be less than 1.";
    //m_data = (T*)malloc(M*N*sizeof(T));
    m_data = new T[M*N];
    for(int i = 0; i < M*N; i++)
        m_data[i] = data[i];
}

template <int M, int N, typename T>
Mat<M,N,T>::~Mat()
{

}

template <int M, int N, typename T>
T* Mat<M,N,T>::data()
{
    return m_data;
}

template <int M, int N, typename T>
T Mat<M,N,T>::at(int i, int j)
{
    return (m_data[i*N+j]);
}

template <int M, int N, typename T>
void Mat<M,N,T>::set(int i, int j, T v)
{
    m_data[i*N+j] = v;
}

template <int M, int N, typename T>
T Mat<M,N,T>::determinant()
{
    if(M!=N)
        return 0;
    return determinantHelper(m_data, M);
}

template <int M, int N, typename T>
T Mat<M,N,T>::determinantHelper(T *data, int dim)
{
    if(dim == 1)
        return data[0];
    if(dim == 2)
        return data[0]*data[3]-data[1]*data[2];
    T deter = 0;
    for(int i = 0; i < dim; i++)
    {
        T* minData = minorData(data, 0, i, dim);
        deter += data[i]*pow(-1,i)*determinantHelper(minData, dim-1);
    }
    return deter;
}

template <int M, int N, typename T>
Mat<N,M,T> Mat<M,N,T>::transposed()
{
    T* newData = (T*)calloc(N*M, sizeof(T));
    for(int i = 0; i < N; i++)
        for(int j = 0; j < M; j++)
            newData[i*M+j] = m_data[j*N+i];
    return (Mat<N,M,T>(newData));
}

template <int M, int N, typename T>
Mat<M,N,T> Mat<M,N,T>::inverse()
{
    if(M!=N)
        throw "Error: M!=N";
    T deter = determinant();
    if(deter == 0)
    {
        cout<<"Error: matrix with determinant=0 has no inverse!\n";
        printData();
        throw "Error: matrix with determinant=0 has no inverse!\n";
    }
    if(M == 1)
    {
        Mat<M,N,T> ret = Mat<M,N,T>();
        ret.set(0, 0, 1.0f/m_data[0]);
        return ret;
    }
    if(M == 2)
    {
        Mat<M,N,T> ret = Mat<M,N,T>();
        ret.set(0, 0, m_data[3]/(deter));
        ret.set(0, 1, m_data[2]/(-deter));
        ret.set(1, 0, m_data[1]/(-deter));
        ret.set(1, 1, m_data[0]/(deter));
        return ret;
    }
    T* data = (T*)malloc(M*N*sizeof(T));
    for(int i = 0; i < M; i++)
        for(int j = 0; j < N; j++)
            data[i*N+j] = pow(-1,i+j)*(cofactor(i, j).determinant())/(deter);
    return Mat<M,N,T>(data).transposed();
}

template <int M, int N, typename T>
Mat<M-1,N-1,T> Mat<M,N,T>::cofactor(int i, int j)
{
    T* data = minorData(m_data, i, j, M);
    Mat<M-1,N-1,T> ret = Mat<M-1,N-1,T>(data);
    return ret;
}


template <int M, int N, typename T>
T* Mat<M,N,T>::minorData(T* olddata, int m, int n, int dim)
{
    T* data = (T*)calloc((dim-1)*(dim-1), sizeof(T));
    int r = 0;
    int c = 0;
    for(int i = 0; i < dim; i++)
    {
        if(i == m) continue;
        c = 0;
        for(int j = 0; j < dim; j++)
        {
            if(j == n) continue;
            data[r*(dim-1)+c] = olddata[i*(dim)+j];
            c++;
        }
        r++;
    }
    return data;
}


#endif // MAT_H
