// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Sum(IList{double})"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Sum", typeof(double[]))]
	[Quality(QualityBand.Mature)]
	public static class FastSumOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'sum'.</param>
		/// <param name="array">Constant value for 'array'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sum,array))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double sum, IList<double> array)
		{
			return (sum == Factor.Sum(array)) ? 0.0 : Double.NegativeInfinity;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'sum'.</param>
		/// <param name="array">Constant value for 'array'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sum,array))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double sum, IList<double> array) { return LogAverageFactor(sum, array); }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sum">Constant value for 'sum'.</param>
		/// <param name="array">Constant value for 'array'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sum,array))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(double sum, IList<double> array) { return LogAverageFactor(sum, array); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Incoming message from 'sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_sum">Outgoing message to 'sum'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sum) p(sum) factor(sum,array))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sum"/> is not a proper distribution</exception>
		public static double LogAverageFactor([SkipIfUniform] Gaussian sum, [Fresh] Gaussian to_sum)
		{
			return to_sum.GetLogAverageOf(sum);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'sum'.</param>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array) p(array) factor(sum,array))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static double LogAverageFactor(double sum, [SkipIfAnyUniform] IList<Gaussian> array)
		{
			Gaussian to_sum = SumAverageConditional(array);
			return to_sum.GetLogProb(sum);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Incoming message from 'sum'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sum) p(sum) factor(sum,array) / sum_sum p(sum) messageTo(sum))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gaussian sum) { return 0.0; }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'sum'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array) p(array) factor(sum,array))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double sum, IList<Gaussian> array) { return LogAverageFactor(sum, array); }

		/// <summary>
		/// EP message to 'sum'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'sum' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sum' as the random arguments are varied.
		/// The formula is <c>proj[p(sum) sum_(array) p(array) factor(sum,array)]/p(sum)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static Gaussian SumAverageConditional([SkipIfAnyUniform] IList<Gaussian> array)
		{
			double mean = 0;
			double variance = 0;
			for (int i = 0; i < array.Count; i++) {
				if (array[i].Precision == 0) return array[i];
				double mean1;
				double variance1;
				array[i].GetMeanAndVarianceImproper(out mean1, out variance1);
				mean = mean + mean1;
				variance = variance + variance1;
			}
			return new Gaussian(mean, variance);
		}

#if true
		/// <summary>
		/// EP message to 'array'
		/// </summary>
		/// <param name="sum">Incoming message from 'sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_sum">Outgoing message to 'sum'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'array' as the random arguments are varied.
		/// The formula is <c>proj[p(array) sum_(sum) p(sum) factor(sum,array)]/p(array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sum"/> is not a proper distribution</exception>
		public static GaussianList ArrayAverageConditional<GaussianList>([SkipIfUniform] Gaussian sum, [Fresh] Gaussian to_sum, IList<Gaussian> array, GaussianList result)
			where GaussianList : IList<Gaussian>, SettableToUniform
		{
			// TM: It is tempting to put SkipIfAllUniform on array but this isn't correct if the array has one element.
			if (array.Count == 1) {
				result[0] = sum;
				return result;
			}
			if (sum.Precision == 0) {
				for (int i = 0; i < result.Count; i++) {
					result[i] = sum;
				}
				return result;
			}
			double totalMean, totalVariance;
			double sumMean, sumVariance;
			sum.GetMeanAndVarianceImproper(out sumMean, out sumVariance);

			// Check if an element of the array is uniform.
			int indexOfUniform = -1;
			for (int i = 0; i < array.Count; i++) {
				double meani, variancei;
				array[i].GetMeanAndVarianceImproper(out meani, out variancei);
				// instead of testing IsUniform, we need to test the more strict requirement 
				// of variance < Infinity due to the way we are doing the computations
				if (double.IsPositiveInfinity(variancei)) {
					if (indexOfUniform >= 0) {
						// more than one element of array is uniform.
						result.SetToUniform();
						return result;
					}
					indexOfUniform = i;
				}
			}
			if (indexOfUniform >= 0) {
				// exactly one element of array is uniform.
				totalMean = 0;
				totalVariance = 0;
				for (int i = 0; i < array.Count; i++) {
					if (i == indexOfUniform) continue;
					double meani, variancei;
					array[i].GetMeanAndVarianceImproper(out meani, out variancei);
					totalMean += meani;
					totalVariance += variancei;
					result[i] = new Gaussian();
				}
				// totalMean = sum_{i except indexOfUniform} array[i].GetMean()
				// totalVariance = sum_{i except indexOfUniform} array[i].GetVariance()
				result[indexOfUniform] = new Gaussian(sumMean - totalMean, sumVariance + totalVariance);
				return result;
			}
			// at this point, the array has no uniform elements.

			// get the mean and variance of sum of all the Gaussians;
			to_sum.GetMeanAndVarianceImproper(out totalMean, out totalVariance);

			// subtract it off from the mean and variance of incoming Gaussian from Sum
			totalMean = sumMean - totalMean;
			totalVariance = sumVariance + totalVariance;

			for (int i = 0; i < array.Count; i++) {
				double meani, variancei;
				array[i].GetMeanAndVarianceImproper(out meani, out variancei);
				result[i] = new Gaussian(totalMean + meani, totalVariance - variancei);
			}
			return result;
		}

		/// <summary>
		/// EP message to 'array'
		/// </summary>
		/// <param name="sum">Incoming message from 'sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_sum">Outgoing message to 'sum'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'array' as the random arguments are varied.
		/// The formula is <c>proj[p(array) sum_(sum) p(sum) factor(sum,array)]/p(array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sum"/> is not a proper distribution</exception>
		public static Gaussian[] ArrayAverageConditional([SkipIfUniform] Gaussian sum, [Fresh] Gaussian to_sum, Gaussian[] array, Gaussian[] result)
		{
			// TM: It is tempting to put SkipIfAllUniform on array but this isn't correct if the array has one element.
			if (array.Length == 1) {
				result[0] = sum;
				return result;
			}
			if (sum.Precision == 0) {
				for (int i = 0; i < result.Length; i++) {
					result[i] = sum;
				}
				return result;
			}
			double totalMean, totalVariance;
			double sumMean, sumVariance;
			sum.GetMeanAndVarianceImproper(out sumMean, out sumVariance);

			// Check if an element of the array is uniform.
			int indexOfUniform = -1;
			for (int i = 0; i < array.Length; i++) {
				double meani, variancei;
				array[i].GetMeanAndVarianceImproper(out meani, out variancei);
				// instead of testing IsUniform, we need to test the more strict requirement 
				// of variance < Infinity due to the way we are doing the computations
				if (double.IsPositiveInfinity(variancei)) {
					if (indexOfUniform >= 0) {
						// more than one element of array is uniform.
						for (int j = 0; j < result.Length; j++) {
							result[j].SetToUniform();							
						}
						return result;
					}
					indexOfUniform = i;
				}
			}
			if (indexOfUniform >= 0) {
				// exactly one element of array is uniform.
				totalMean = 0;
				totalVariance = 0;
				for (int i = 0; i < array.Length; i++) {
					if (i == indexOfUniform) continue;
					double meani, variancei;
					array[i].GetMeanAndVarianceImproper(out meani, out variancei);
					totalMean += meani;
					totalVariance += variancei;
					result[i] = new Gaussian();
				}
				// totalMean = sum_{i except indexOfUniform} array[i].GetMean()
				// totalVariance = sum_{i except indexOfUniform} array[i].GetVariance()
				result[indexOfUniform] = new Gaussian(sumMean - totalMean, sumVariance + totalVariance);
				return result;
			}
			// at this point, the array has no uniform elements.

			// get the mean and variance of sum of all the Gaussians;
			to_sum.GetMeanAndVarianceImproper(out totalMean, out totalVariance);

			// subtract it off from the mean and variance of incoming Gaussian from Sum
			totalMean = sumMean - totalMean;
			totalVariance = sumVariance + totalVariance;

			for (int i = 0; i < array.Length; i++) {
				double meani, variancei;
				array[i].GetMeanAndVarianceImproper(out meani, out variancei);
				Gaussian item = new Gaussian(totalMean + meani, totalVariance - variancei);
				result[i] = item;
			}
			return result;
		}

		/// <summary>
		/// EP message to 'array'
		/// </summary>
		/// <param name="sum">Constant value for 'sum'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' conditioned on the given values.
		/// </para></remarks>
		public static GaussianList ArrayAverageConditional<GaussianList>([SkipIfUniform] double sum, IList<Gaussian> array, GaussianList result)
			where GaussianList : IList<Gaussian>, SettableToUniform
		{
			// TM: It is tempting to put SkipIfAllUniform on array but this isn't correct if the array has one element.
			Gaussian to_sum = SumAverageConditional(array);
			return ArrayAverageConditional(Gaussian.PointMass(sum), to_sum, array, result);
		}

#else
		public static Gaussian ArrayAverageConditional([SkipIfUniform] Gaussian sum, [Fresh] Gaussian to_sum, [Proper, MatchingIndex] IList<Gaussian> array, int resultIndex)
		{
		// TM: It is tempting to put SkipIfAllUniform on array but this isn't correct if the array has one element.
			if (sum.Precision == 0) return sum;
			double sumMean, sumVar;
			sum.GetMeanAndVarianceImproper(out sumMean, out sumVar);

			double arraySumOfMean, arraySumOfVar;
			if (array[resultIndex].Precision == 0) {
				arraySumOfMean = 0;
				arraySumOfVar = 0;
				for (int i = 0; i < array.Count; i++) {
					if (i == resultIndex) continue;
					if (array[i].Precision == 0) return array[i];
					double meani, variancei;
					array[i].GetMeanAndVarianceImproper(out meani, out variancei);
					arraySumOfMean += meani;
					arraySumOfVar += variancei;
				}
			} else {
				to_sum.GetMeanAndVariance(out arraySumOfMean, out arraySumOfVar);
				double meani, variancei;
				array[resultIndex].GetMeanAndVarianceImproper(out meani, out variancei);
				arraySumOfMean -= meani;
				arraySumOfVar -= variancei;
			}
			return new Gaussian(sumMean - arraySumOfMean, sumVar + arraySumOfVar);
		}

		public static Gaussian ArrayAverageConditional([SkipIfUniform] double sum, [Proper,MatchingIndex] IList<Gaussian> array, int resultIndex)
		{
			Gaussian to_sum = SumAverageConditional(array);
			return ArrayAverageConditional(Gaussian.PointMass(sum), to_sum, array, resultIndex);
		}
#endif

		// VMP //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/// <summary>
		/// VMP message to 'sum'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'sum' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sum' as the random arguments are varied.
		/// The formula is <c>proj[sum_(array) p(array) factor(sum,array)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static Gaussian SumAverageLogarithm([SkipIfAnyUniform] IList<Gaussian> array)
		{
			return SumAverageConditional(array);
		}

		public static Gaussian SumDeriv(Gaussian sum, [Proper] IList<Gaussian> array, IList<Gaussian> to_array)
		{
			double sumPrec = sum.Precision;
			double sumOfArrayVariance = 0;
			double sumOfPostVariance = 0;
			for (int i = 0; i < array.Count; i++) {
				double mi,vi;
				Gaussian prior = array[i]/to_array[i];
				prior.GetMeanAndVariance(out mi, out vi);
				sumOfArrayVariance += vi;
				double vpost = 1/(sumPrec + prior.Precision);
				sumOfPostVariance += vpost;
			}
			double mpDeriv = 1/(sumPrec + 1/sumOfArrayVariance)/sumOfPostVariance;
			double precDeriv = 1;
			return Gaussian.FromNatural(mpDeriv-1, precDeriv-1);
		}

#if true
		public static GaussianList ArrayAverageLogarithm6<GaussianList>([SkipIfUniform] Gaussian sum, [Proper] IList<Gaussian> array, IList<Gaussian> array_deriv, GaussianList to_array)
			where GaussianList : IList<Gaussian>
		{
			GaussianList result = to_array;
			double sumMean, sumVar;
			sum.GetMeanAndVariance(out sumMean, out sumVar);
			if (to_array[0].IsPointMass) return result;
			if (to_array[0].IsUniform()) return ArrayAverageLogarithm5(sum, array, to_array);
			//if (array_deriv[0].IsUniform()) return ArrayAverageLogarithm2(sum, array, to_array);

			int n = array.Count;
			if (false) { // debugging
				Console.Write("mp is");
				for (int i = 0; i < n; i++) {
					Console.Write(" "+array[i].MeanTimesPrecision);
				}
				Console.WriteLine();
				Console.Write("new mp should be");
				for (int i = 0; i < n; i++) {
					double deriv = array_deriv[i].MeanTimesPrecision+1;
					Console.Write(" "+(array[i].MeanTimesPrecision+deriv));
					double tau = to_array[i].MeanTimesPrecision + 1;
					result[i] = new Gaussian(tau*sumVar, sumVar);
				}
				Console.WriteLine();
				return result;
			}
			double bsum = 0, asum = 0;
			for (int i = 0; i < array.Count; i++) {
				double mi,vi;
				array[i].GetMeanAndVariance(out mi, out vi);
				double deriv = array_deriv[i].MeanTimesPrecision+1;
				double offset = array[i].MeanTimesPrecision - deriv*to_array[i].MeanTimesPrecision;
				//prior.Precision = array[i].Precision - array_deriv[i].Precision*to_array[i].Precision;
				double b = vi*(offset + deriv*sum.MeanTimesPrecision);
				double a = vi*deriv*sum.Precision;
				double r = 1/(1-a);
				bsum += b*r;
				asum += a*r;
			}
			double c = bsum/(1 + asum);
			for (int i = 0; i < array.Count; i++) {
				double mi,vi;
				array[i].GetMeanAndVariance(out mi, out vi);
				double deriv = array_deriv[i].MeanTimesPrecision+1;
				double offset = array[i].MeanTimesPrecision - deriv*to_array[i].MeanTimesPrecision;
				double b = vi*(offset + deriv*sum.MeanTimesPrecision);
				double a = vi*deriv*sum.Precision;
				double mnew = (b-a*c)/(1-a);
				double tau = (mnew/vi - offset)/deriv;
				result[i] = new Gaussian(tau*sumVar, sumVar);
			}
			return result;
		}
		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="sum">Incoming message from 'sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="to_array">Previous outgoing message to 'array'.</param>
		/// <returns>The outgoing VMP message to the 'array' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' with 'sum' integrated out.
		/// The formula is <c>sum_sum p(sum) factor(sum,array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sum"/> is not a proper distribution</exception>
		public static GaussianList ArrayAverageLogarithm2<GaussianList>([SkipIfUniform] Gaussian sum, [Proper] IList<Gaussian> array, GaussianList to_array)
			where GaussianList : IList<Gaussian>
		{
			GaussianList result = to_array;
			double sumMean, sumVar;
			sum.GetMeanAndVariance(out sumMean, out sumVar);
			if (to_array[0].IsPointMass) return result;

			// Since the approximate marginals q(array[i]) are coupled, instead of updating them once, solve for their converged value directly.
			double sumOfArrayMean = 0;
			double sumOfArrayVariance = 0;
			for (int i = 0; i < array.Count; i++) {
				double mi,vi;
				Gaussian prior = array[i]/to_array[i];
				prior.GetMeanAndVariance(out mi, out vi);
				sumOfArrayMean += mi;
				sumOfArrayVariance += vi;
			}
			double totalVariance = sumVar + sumOfArrayVariance;
			if (totalVariance == 0) return result;
			double r = (sumMean - sumOfArrayMean) / totalVariance;
			for (int i = 0; i < result.Count; i++) {
				double mi,vi;
				Gaussian prior = array[i]/to_array[i];
				prior.GetMeanAndVariance(out mi, out vi);
				double mean = r*(sumVar + vi) + mi;
				result[i] = new Gaussian(mean, sumVar);
			}
			return result;
		}
		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="sum">Incoming message from 'sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="to_array">Previous outgoing message to 'array'.</param>
		/// <returns>The outgoing VMP message to the 'array' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' with 'sum' integrated out.
		/// The formula is <c>sum_sum p(sum) factor(sum,array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sum"/> is not a proper distribution</exception>
		public static GaussianList ArrayAverageLogarithm1<GaussianList>([SkipIfUniform] Gaussian sum, [Stochastic, Proper] IList<Gaussian> array, GaussianList to_array)
			where GaussianList : IList<Gaussian>
		{
			GaussianList result = to_array;
			double sumMean, sumVar;
			sum.GetMeanAndVariance(out sumMean, out sumVar);

			// This version does one update of q(array[i]) for each array element in turn.
			double arraySumOfMean = 0;
			for (int i = 0; i < array.Count; i++) {
				arraySumOfMean = arraySumOfMean + array[i].GetMean();
			}
			double stepsize = 1;
			for (int i = 0; i < result.Count; i++) {
				double partialSum = arraySumOfMean - array[i].GetMean();
				Gaussian oldresult = result[i];
				result[i] = new Gaussian(sumMean - partialSum, sumVar);
				Gaussian newMarginal = array[i] * ((result[i] / oldresult) ^ stepsize);
				arraySumOfMean = partialSum + newMarginal.GetMean();
			}
			return result;
		}
		public static GaussianList ArrayAverageLogarithm5<GaussianList>([SkipIfUniform] Gaussian sum, [Proper] IList<Gaussian> array, GaussianList to_array)
			where GaussianList : IList<Gaussian>
		{
			GaussianList result = to_array;
			double sumMean, sumVar;
			sum.GetMeanAndVariance(out sumMean, out sumVar);

			// This version does one update of q(array[i]) for each array element in turn.
			double arraySumOfMean = 0;
			for (int i = 0; i < array.Count; i++) {
				arraySumOfMean = arraySumOfMean + array[i].GetMean();
			}
			for (int i = 0; i < result.Count; i++) {
				double partialSum = arraySumOfMean - array[i].GetMean();
				Gaussian oldresult = result[i];
				result[i] = new Gaussian(sumMean - partialSum, sumVar);
			}
			return result;
		}
		public static GaussianList ArrayAverageLogarithm4<GaussianList>([SkipIfUniform] Gaussian sum, [Proper] IList<Gaussian> array, GaussianList to_array)
			where GaussianList : IList<Gaussian>, ICloneable
		{
			GaussianList array2 = (GaussianList)to_array.Clone();
			for (int i = 0; i < array.Count; i++) {
				array2[i] = array[i];
			}
			GaussianList result = (GaussianList)to_array.Clone();
			for (int iter = 0; iter < 10; iter++) {
				GaussianList oldresult = (GaussianList)result.Clone();
				result = ArrayAverageLogarithm1(sum, array2, result);
				for (int i = 0; i < array2.Count; i++) {
					array2[i] = array2[i]*(result[i]/oldresult[i]);
				}
			}
			return result;
		}
		public static GaussianList ArrayAverageLogarithm3<GaussianList>([SkipIfUniform] Gaussian sum, [Proper] IList<Gaussian> array, GaussianList to_array)
			where GaussianList : IList<Gaussian>, ICloneable
		{
			//GaussianList result1 = ArrayAverageLogarithm(sum, array, to_array);
			//ArrayAverageLogarithm1(sum, array, to_array);
			GaussianList result = to_array;
			double sumMean, sumVar;
			sum.GetMeanAndVariance(out sumMean, out sumVar);

			double arraySumOfMean = 0;
			for (int i = 0; i < array.Count; i++) {
				arraySumOfMean = arraySumOfMean + array[i].GetMean();
			}

			int n = array.Count;
			Vector a = Vector.Zero(n);
			Vector b = Vector.Zero(n);
			Gaussian[] newMarginal = new Gaussian[n];
			for (int i = 0; i < n; i++) {
				double mi,vi;
				array[i].GetMeanAndVariance(out mi, out vi);
				double partialSum = arraySumOfMean - mi;
				double mr,vr;
				to_array[i].GetMeanAndVariance(out mr, out vr);
				double mu1 = sumMean - partialSum;
				newMarginal[i] = array[i]*(new Gaussian(mu1, sumVar)/to_array[i]);
				double vpost = 1/(1/vi + 1/sumVar - 1/vr);
				double m1 = vpost*(mi/vi + mu1/sumVar - mr/vr);
				arraySumOfMean = partialSum + m1;
				a[i] = -vpost/sumVar;
				b[i] = m1 + vpost*(sumMean/sumVar) - vpost*mu1/sumVar;
			}
			Matrix x = new Matrix(n, n);
			for (int i = 0; i < n; i++) {
				x[i, i] = 1 + a[i];
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					x[i, j] -= a[i];
				}
			}
			(new LuDecomposition(x)).Solve(b);
			arraySumOfMean = 0;
			for (int i = 0; i < n; i++) {
				arraySumOfMean += b[i];
			}
			for (int i = 0; i < n; i++) {
				result[i] = new Gaussian(sumMean - arraySumOfMean + b[i], sumVar);
			}
			return result;
		}

		public static GaussianList ArrayAverageLogarithm<GaussianList>([SkipIfUniform] Gaussian sum, [Proper] IList<Gaussian> array, GaussianList to_array)
			where GaussianList : IList<Gaussian>, ICloneable
		{
			return ArrayAverageLogarithm1(sum, array, to_array);
		}
		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="sum">Constant value for 'sum'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' conditioned on the given values.
		/// </para></remarks>
		public static GaussianList ArrayAverageLogarithm<GaussianList>(double sum, [Proper] IList<Gaussian> array, GaussianList result)
			where GaussianList : IList<Gaussian>, ICloneable
		{
			return ArrayAverageLogarithm(Gaussian.PointMass(sum), array, result);
		}
#else
		public static Gaussian ArrayAverageLogarithm([SkipIfUniform] Gaussian sum, [Fresh] Gaussian to_sum, [Proper, MatchingIndex] IList<Gaussian> array, int resultIndex)
		{
			double sumMean, sumVar;
			sum.GetMeanAndVariance(out sumMean, out sumVar);

			double arraySumOfMean, arraySumOfVar;
			to_sum.GetMeanAndVariance(out arraySumOfMean, out arraySumOfVar);
			arraySumOfMean -= array[resultIndex].GetMean();
			return new Gaussian(sumMean - arraySumOfMean, sumVar);
		}
#endif

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sum,array))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }
	}

	/////////////////////////////////////////////////////////////////////////////
	// Example Code for user guide
	//[FactorMethod(typeof(Factor), "Sum", typeof(double[]))]
#if SUPPRESS_XMLDOC_WARNINGS
#pragma warning disable 1591
#endif
	internal static class SumOp
	{
		public static Gaussian SumAverageConditional([SkipIfAnyUniform] IList<Gaussian> array)
		{
			double mean=0;
			double variance=0;
			double mean1;
			double variance1;
			for (int i = 0; i < array.Count; i++) {
				array[i].GetMeanAndVariance(out mean1, out variance1);
				mean = mean + mean1;
				variance = variance + variance1;
			}
			return new Gaussian(mean, variance);
		}

		public static GaussianArray ArrayAverageConditional<GaussianArray>([SkipIfAnyUniform] GaussianArray array, [SkipIfUniform] Gaussian sum, GaussianArray result)
			where GaussianArray : IList<Gaussian>
		{
			double mean, mean1;
			double variance, variance1;
			// get the mean and variance of sum of all the Gaussians
			Gaussian to_sum = SumAverageConditional(array);
			to_sum.GetMeanAndVariance(out mean, out variance);
			// subtract it off from the mean and variance of incoming Gaussian from Sum
			sum.GetMeanAndVariance(out mean1, out variance1);
			mean = mean1 - mean;
			variance = variance1 + variance;
			for (int i = 0; i < array.Count; i++) {
				array[i].GetMeanAndVariance(out mean1, out variance1);
				result[i] = new Gaussian(mean + mean1, variance - variance1);
			}
			return result;
		}

		public static GaussianArray ArrayAverageConditional<GaussianArray>([SkipIfAnyUniform] GaussianArray array, double sum, GaussianArray result)
			where GaussianArray : IList<Gaussian>
		{
			return ArrayAverageConditional(array, Gaussian.PointMass(sum), result);
		}

		public static double LogAverageFactor(double sum, [SkipIfAnyUniform] IList<Gaussian> array)
		{
			Gaussian to_sum = SumAverageConditional(array);
			return to_sum.GetLogProb(sum);
		}
		public static double LogEvidenceRatio(double sum, IList<Gaussian> array) { return LogAverageFactor(sum, array); }

		public static double LogAverageFactor([SkipIfUniform] Gaussian sum, [SkipIfAnyUniform] IList<Gaussian> array)
		{
			Gaussian to_sum = SumAverageConditional(array);
			return to_sum.GetLogAverageOf(sum);
		}
		[Skip]
		public static double LogEvidenceRatio(Gaussian sum) { return 0.0; }

	}
}


