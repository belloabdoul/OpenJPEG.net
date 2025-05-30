﻿using System;

namespace OpenJpeg.Internal;

internal static class Invert
{
	/**
	 * Matrix inversion.
	 */
	public static bool MatrixInversion(float[] SrcMatrix,
		float[] DestMatrix,
		int n_compo)
	{
		var permutation_size = n_compo;
		var swap_size = n_compo;
            
		var permutations = GC.AllocateUninitializedArray<int>(permutation_size);
		var double_data = GC.AllocateUninitializedArray<float>(swap_size);
		var dest_temp = GC.AllocateUninitializedArray<float>(swap_size);
		var swap_area = GC.AllocateUninitializedArray<float>(swap_size);

		if (!lupDecompose(SrcMatrix, permutations, double_data, n_compo))
			return false;

		lupInvert(SrcMatrix, DestMatrix, n_compo, permutations, double_data, dest_temp, swap_area);

		return true;
	}

	private static bool lupDecompose(float[] matrix, int[] permutations, 
		float[] swap_area, int n_compo) 
	{
		var tmpPermutations = 0;
		var k2=0;
		float temp;
		int i;
		var lLastColum = n_compo - 1;
		var lSwapSize = n_compo;
		var lTmpMatrix = 0;
		var offset = 1;
		var lStride = n_compo-1;

		/*initialize permutations */
		for (i = 0; i < n_compo; ++i) 
			permutations[tmpPermutations++] = i;
	        
		/* now make a pivot with colum switch */
		tmpPermutations = 0;
		for (var k = 0; k < lLastColum; ++k) 
		{
			var p = 0f;

			/* take the middle element */
			var lColumnMatrix = lTmpMatrix + k;

			/* make permutation with the biggest value in the column */
			for (i = k; i < n_compo; ++i) 
			{
				temp = matrix[lColumnMatrix] > 0 ? matrix[lColumnMatrix] : -matrix[lColumnMatrix];
				if (temp > p) {
					p = temp;
					k2 = i;
				}

				/* next line */
				lColumnMatrix += n_compo;
			}

			/* a whole rest of 0 -> non singular */
			if (p == 0.0)
				return false;

			/* should we permute ? */
			if (k2 != k) {
				/*exchange of line */
				/* k2 > k */
				var dstPermutations = tmpPermutations + k2 - k;
				/* swap indices */
				var t = permutations[tmpPermutations];
				permutations[tmpPermutations] = permutations[dstPermutations];
				permutations[dstPermutations] = t;

				/* and swap entire line. */
				lColumnMatrix = lTmpMatrix + (k2 - k) * n_compo;
				Array.Copy(matrix, lColumnMatrix, swap_area, 0, n_compo);
				Array.Copy(matrix, lTmpMatrix, matrix, lColumnMatrix, n_compo);
				Array.Copy(swap_area, 0, matrix, lTmpMatrix, n_compo);
			}

			/* now update data in the rest of the line and line after */
			var lDestMatrix = lTmpMatrix + k;
			lColumnMatrix = lDestMatrix + n_compo;
			/* take the middle element */
			temp = matrix[lDestMatrix++];

			/* now compute up data (i.e. coeff up of the diagonal). */
			for (i = offset; i < n_compo; ++i)  {
				/*lColumnMatrix; */
				/* divide the lower column elements by the diagonal value */
				/* matrix[i][k] /= matrix[k][k]; */
				/* p = matrix[i][k] */
				p = matrix[lColumnMatrix] / temp;
				matrix[lColumnMatrix++] = p;
     		
				for ( /* k + 1 */var j = offset; j < n_compo; ++j) {
					/* matrix[i][j] -= matrix[i][k] * matrix[k][j]; */
					matrix[lColumnMatrix++] -= p * matrix[lDestMatrix++];
				}

				/* come back to the k+1th element */
				lDestMatrix -= lStride;
				/* go to kth element of the next line */
				lColumnMatrix += k;
			}

			/* offset is now k+2 */
			++offset;
			/* 1 element less for stride */
			--lStride;
			/* next line */
			lTmpMatrix+=n_compo;
			/* next permutation element */
			++tmpPermutations;
		}
		return true;
	}

	private static void lupSolve(float[] result,
		float[] matrix,
		float[] vector,
		int[] permutations,
		int n_compo, float[] intermediate_data)
	{
		int j;
		float sum;
		var lStride = n_compo + 1;
		int lCurrentPtr;
		int lTmpMatrix;
		var lLineMatrix = 0;
		var lBeginPtr = n_compo - 1;
		var lCurrentPermutationPtr = 0;


		var lIntermediatePtr = 0;
		var lGeneratedData = n_compo - 1;

		for (var i = 0; i < n_compo; ++i)
		{
			sum = 0f;
			lCurrentPtr = 0;
			lTmpMatrix = lLineMatrix;
			for (j = 1; j <= i; ++j)
			{
				/* sum += matrix[i][j-1] * y[j-1]; */
				sum += matrix[lTmpMatrix++] * intermediate_data[lCurrentPtr++];
			}

			/*y[i] = pVector[pPermutations[i]] - sum; */
			intermediate_data[lIntermediatePtr++] = vector[permutations[lCurrentPermutationPtr++]] - sum;
			lLineMatrix += n_compo;
		}

		/* we take the last point of the matrix */
		lLineMatrix = n_compo * n_compo - 1;

		/* and we take after the last point of the destination vector */
		var lDestPtr = n_compo;


		System.Diagnostics.Debug.Assert(n_compo != 0);
		for (var k = n_compo - 1; k != -1; --k)
		{
			sum = 0f;
			lTmpMatrix = lLineMatrix;
			var u = matrix[lTmpMatrix++];
			lCurrentPtr = lDestPtr--;
			for (j = k + 1; j < n_compo; ++j)
			{
				/* sum += matrix[k][j] * x[j] */
				sum += matrix[lTmpMatrix++] * intermediate_data[lCurrentPtr++];
			}

			/*x[k] = (y[k] - sum) / u; */
			result[lBeginPtr--] = (intermediate_data[lGeneratedData--] - sum) / u;
			lLineMatrix -= lStride;
		}
	}

	private static void lupInvert(float[] SrcMatrix,
		float[] DestMatrix,
		int n_compo,
		int[] permutations,
		float[] src_temp,
		float[] dest_temp,
		float[] swap_area)
	{
		var lLineMatrix = 0;

		for (var j = 0; j < n_compo; ++j)
		{
			var lCurrentPtr = lLineMatrix++;
			for (var c = 0; c < src_temp.Length; c++)
				src_temp[c] = 0;
			src_temp[j] = 1f;
			lupSolve(dest_temp, SrcMatrix, src_temp, permutations, n_compo, swap_area);

			for (var i = 0; i < n_compo; ++i)
			{
				DestMatrix[lCurrentPtr] = dest_temp[i];
				lCurrentPtr += n_compo;
			}
		}
	}
}