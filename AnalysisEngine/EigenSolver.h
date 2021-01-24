//#pragma once
//#include "eigen/sparsecholesky"
//#include "eigen"
//#include "eigen/sparse"
//#include <vector>
//
//class eigensolver
//{
//public:
//	typedef eigen::sparsematrix<double> spmat; // declares a column-major sparse matrix type of double
//	typedef eigen::triplet<double> t;
//
//	static std::vector<double> solve(std::vector<std::vector<double>> a, std::vector<double> b)
//	{
//	
//		int n = a.size();
//
//		// assembly:
//		eigen::vectorxd bbar(n);                   // the right hand side-vector resulting from the constraints
//
//		eigen::sparsematrix<double> asd(n, n);
//
//		//for (size_t i = 0; i < n; i++)
//		//{
//		//	asd.push
//		//}
//
//		//// solving:
//		//eigen::simplicialcholesky<spmat> chol(asd);  // performs a cholesky factorization of a
//		//eigen::vectorxd x = chol.solve(bbar);         // use the factorization to solve for the given right hand side
//
//
//		//return null;
//		////return retval;
//	};
//
//private:
//};
//
//
