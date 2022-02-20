#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <random>
using Matrix = std::vector<std::vector<double>>;
namespace algebra{
Matrix zeros(size_t n, size_t m){
    if (n<=0 || m<=0)
        throw std::logic_error("n and m can not be less than 0");
    Matrix zero_mat(n);
    for(size_t i{};i<n;i++)
        zero_mat[i].resize(m);
    return zero_mat;
}
Matrix ones(size_t n, size_t m){
    if (n<=0 || m<=0)
        throw std::logic_error("n and m can not be less than 0");
    Matrix ones_mat(n);
    for(size_t i{};i<n;i++)
        ones_mat[i].resize(m);
    for(size_t i{};i<n;i++)
        for(size_t j{};j<m;j++)
        ones_mat[i][j]=1;
    return ones_mat;
}
Matrix random(size_t n, size_t m, double min, double max){
    if (n<=0 || m<=0)
        throw std::logic_error("n and m can not be less than 0");
    if (min>max)
        throw std::logic_error("min can not be greater than min");
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(min, max);
    Matrix rand_mat=zeros(n,m);
    for(size_t i{};i<n;i++)
        for(size_t j{};j<m;j++)
            rand_mat[i][j]=dist(mt);
    return rand_mat;
}
void show(const Matrix& matrix){
    for(size_t i{};i<matrix.size();i++){
        for(size_t j{};j<matrix[0].size();j++)
            std::cout<<std::setprecision(3)<<std::setw(10)<<matrix[i][j];
        std::cout<<std::endl;}
}
Matrix multiply(const Matrix& matrix, double c){
    Matrix multi_mat=zeros(matrix.size(),matrix[0].size());
    for(size_t i{};i<matrix.size();i++)
        for(size_t j{};j<matrix[0].size();j++)
            multi_mat[i][j]=c*matrix[i][j];
    return multi_mat;
}
Matrix multiply(const Matrix& matrix1, const Matrix& matrix2){
    if (matrix1.empty() || matrix2.empty())
        return{};
    Matrix m=zeros(matrix1.size(),matrix2[0].size());
    if (matrix1[0].size()!=matrix2.size())
        throw std::logic_error("can not be multiplied");
    
    for(size_t i{};i<matrix1.size();i++)
        for(size_t j{};j<matrix2[0].size();j++)
            for(size_t k{};k<matrix2.size();k++)
                m[i][j]+=matrix1[i][k]*matrix2[k][j];
        
        
    return m;
}
Matrix sum(const Matrix& matrix, double c){
    if (!matrix.empty()){
        Matrix scal_sum=zeros(matrix.size(),matrix[0].size());
        for(size_t i{};i<matrix.size();i++)
            for(size_t j{};j<matrix[0].size();j++)
                scal_sum[i][j]=c+matrix[i][j];
        return scal_sum;}
    else
        return matrix;

}
Matrix sum(const Matrix& matrix1, const Matrix& matrix2){
    if (!matrix1.empty()&&!matrix2.empty()){
        if((matrix1.size()==matrix2.size()) && (matrix1[0].size()==matrix2[0].size())){
            Matrix m=zeros(matrix1.size(),matrix1[0].size());
                for(size_t i{};i<matrix1.size();i++)
                    for(size_t j{};j<matrix1[0].size();j++)
                        m[i][j]=matrix1[i][j]+matrix2[i][j];

            return m;}
        else
            throw std::logic_error("two matrix with different dim");}
    else
        {
        if (!matrix1.empty()||!matrix2.empty())
             throw std::logic_error("two matrix with different dim");
        return matrix2;}


}
Matrix transpose(const Matrix& matrix){
    if (!matrix.empty()){
        Matrix m=zeros(matrix[0].size(),matrix.size());
            for(size_t i{};i<matrix[0].size();i++)
                for(size_t j{};j<matrix.size();j++)
                    m[i][j]=matrix[j][i];
        return m;}
    else 
        return matrix;
}
Matrix minor(const Matrix& matrix, size_t n, size_t m){
    int k{-1};
    int l{-1};
    size_t n1=matrix.size();
    size_t m1=matrix[0].size();
    if (n>n1-1 || m>m1-1)
        throw std::logic_error("input rows or coloumns are greater than matrix size");
    Matrix m2=zeros(n1-1,m1-1);
        for(size_t i{};i<n1;i++){
            if (i==n)
                continue;
            else
                {k++;
            for(size_t j=0;j<m1;j++){
                if(j==m)
                    continue;
                else
                    {l++;
                    m2[k][l]=matrix[i][j];}}
                    l=-1;}}
    return m2;
}
double determinant(const Matrix& matrix){
    if (matrix.empty())
        return 1;
    else if(!(matrix.size()==matrix[0].size()))
        throw std::logic_error("the dims are not same");
    else{
    double det_mat{0.0};
    if (matrix.size()==2)
        det_mat=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0];
    else
        for (size_t i{};i<matrix.size();i++)
            det_mat+=matrix[0][i]*pow(-1,i)*determinant(minor(matrix,0,i));
    return det_mat;}
}
 Matrix inverse(const Matrix& matrix){
    if (matrix.empty())
        return matrix;
    if(!(matrix.size()==matrix[0].size()))
        throw std::logic_error("the dims are not same");
    if (matrix.empty())
        return matrix;
    else if (determinant(matrix)==0)
        throw std::logic_error("singular matrix has not inverse");
    else{
    Matrix inv_mat=zeros(matrix.size(),matrix[0].size());
    for(size_t i{};i<matrix.size();i++){
            for(size_t j{};j<matrix.size();j++)
                inv_mat[i][j]=pow(-1,i+j)*determinant(minor(matrix,i,j));}
    inv_mat=multiply(transpose(inv_mat),1/determinant(matrix));
    return inv_mat;}

}
 Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2,int axis){
    size_t n1=matrix1.size();
    size_t m1=matrix1[0].size();
    size_t n2=matrix2.size();
    size_t m2=matrix2[0].size();
    Matrix matrix1_co;
    if (axis==0){
        if (!(m1==m2))
            throw std::logic_error("wrong dims");
        else
        {
            matrix1_co=zeros(n1+n2,m1);
            for(size_t i=0;i<n1;i++)
                for(size_t j=0;j<m1;j++)
                    matrix1_co[i][j]=matrix1[i][j];

            for(size_t i=n1;i<n1+n2;i++)
                for(size_t j=0;j<m2;j++)
                    matrix1_co[i][j]=matrix2[i-n1][j];
        }
    }
    else{
        if (!(n1==n2))
            throw std::logic_error("wrong dims");
        else
        {  
            matrix1_co=zeros(n1,m1+m2);
            for(size_t i=0;i<n1;i++)
                for(size_t j=0;j<m1;j++)
                    matrix1_co[i][j]=matrix1[i][j];

            for(size_t i=0;i<n2;i++)
                for(size_t j=m1;j<m2+m1;j++)
                    matrix1_co[i][j]=matrix2[i][j-m1];
        }
    }
    return matrix1_co;
 }
 Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2){
    if (matrix.empty())
        return{};
    double temp;
    size_t n=matrix.size();
    size_t m=matrix[0].size();
    Matrix matrix_co=zeros(n,m);
    if ((r1>n-1)||(r2>n-1))
        throw std::logic_error("wrong rows are selected");
    else{
    for(size_t i=0;i<n;i++)
        for(size_t j=0;j<m;j++)
            matrix_co[i][j]=matrix[i][j];
    for (size_t i{};i<m;i++){
        temp=matrix_co[r1][i];
        matrix_co[r1][i]=matrix_co[r2][i];
        matrix_co[r2][i]=temp;}
    return matrix_co;
    }
 }
 Matrix ero_multiply(const Matrix& matrix, size_t r, double c){
    if (matrix.empty())
        throw std::logic_error("input matrix is empty");
    if (r>matrix.size()-1)
        throw std::logic_error("wrong row is selected");
    size_t n=matrix.size();
    size_t m=matrix[0].size();
    Matrix matrix_co=zeros(n,m);
    for(size_t i{};i<n;i++)
        for(size_t j{};j<m;j++)
            matrix_co[i][j]=matrix[i][j];
    for (size_t i{};i<m;i++)
        matrix_co[r][i]=matrix_co[r][i]*c;
    return matrix_co;
 }
Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2){
    if (matrix.empty())
        return{};
    //if (r1>matrix.size()-1 || r2>matrix[0].size()-1)
        //throw std::logic_error("wrong rows are selected");
    size_t n{matrix.size()};
    size_t m{matrix[0].size()};
    Matrix matrix_co{zeros(n,m)};
    for(size_t i{};i<n;i++)
        for(size_t j{};j<m;j++)
            matrix_co[i][j]=matrix[i][j];
    for (size_t i{};i<m;i++)
        matrix_co[r2][i]=matrix_co[r1][i]*c+matrix_co[r2][i];
    return matrix_co;}
Matrix upper_triangular(const Matrix& matrix){
     if (matrix.empty())
        throw std::logic_error("input matrix is empty");
    size_t n{matrix.size()};
    size_t m{matrix[0].size()};
    int count{};
    Matrix matrix_co;
    if (matrix.empty())
        return matrix;
    else if(n!=m)
        throw std::logic_error("matrix should be squre matrix");
    else{
        matrix_co=zeros(n,m);
        for(size_t i=0;i<n;i++)
            for(size_t j=0;j<m;j++)
                matrix_co[i][j]=matrix[i][j];
        for(size_t i=0;i<n;i++){
            int k=i;
            int count1{};
            while(matrix_co[k][count]==0){
                if (matrix_co[k+count1][count]!=0)
                    matrix_co=ero_swap(matrix_co,k,k+count1);
                else
                    count1++;
                }
            for (size_t j=i+1;j<n;j++)
                matrix_co=ero_sum(matrix_co,i,-matrix_co[j][count]/matrix_co[i][count],j);
            count++;}
        return matrix_co;
    }
}
}

