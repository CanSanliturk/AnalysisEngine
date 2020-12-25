//#pragma once
//#include "Eigen/SparseCholesky"
//#include "Eigen"
//#include "Eigen/Sparse"
//#include <vector>
//
//class EigenSolver
//{
//public:
//	typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
//	typedef Eigen::Triplet<double> T;
//
//	static std::vector<double> Solve(std::vector<std::vector<double>> A, std::vector<double> b)
//	{
//	
//		int n = A.size();
//
//		// Assembly:
//		Eigen::VectorXd bbar(n);                   // the right hand side-vector resulting from the constraints
//
//		Eigen::SparseMatrix<double> asd(n, n);
//
//		//for (size_t i = 0; i < n; i++)
//		//{
//		//	asd.push
//		//}
//
//		//// Solving:
//		//Eigen::SimplicialCholesky<SpMat> chol(asd);  // performs a Cholesky factorization of A
//		//Eigen::VectorXd x = chol.solve(bbar);         // use the factorization to solve for the given right hand side
//
//
//		//return null;
//		////return retVal;
//	};
//
//private:
//};
//
//
