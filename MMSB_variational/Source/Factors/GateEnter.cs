// (C) Copyright 2008 Microsoft Research Cambridge
#define UseRatioDir
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;
using MicrosoftResearch.Infer.Collections;

namespace MicrosoftResearch.Infer.Factors
{

	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.EnterPartial{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Gate), "EnterPartial<>", null, typeof(int), null, typeof(int[]))]
	[FactorMethod(typeof(Gate), "EnterPartial<>", null, typeof(bool), null, typeof(int[]))]
	[Quality(QualityBand.Mature)]
	public static class GateEnterPartialOp<T>
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enterPartial,cases,value,indices))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio()
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enterPartial,cases,value,indices))</c>.
		/// </para></remarks>
		[Skip]
		public static double LogAverageFactor()
		{
			return 0.0;
		}

		/// <summary>
		/// EP message to 'enterPartial'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'enterPartial' as the random arguments are varied.
		/// The formula is <c>proj[p(enterPartial) sum_(value) p(value) factor(enterPartial,cases,value,indices)]/p(enterPartial)</c>.
		/// </para></remarks>
		public static TList EnterPartialAverageConditional<T, TList>([IsReturnedInEveryElement] T value, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(value);
			return result;
		}
		[Skip]
		public static ArrayType EnterPartialInit<T, ArrayType>([IgnoreDependency] T value, int[] indices, IArrayFactory<T,ArrayType> factory)
			where T : ICloneable
		{
			return factory.CreateArray(indices.Length, i => (T)value.Clone());
		}

		[Skip]
		public static Discrete SelectorAverageConditional(Discrete result)
		{
			result.SetToUniform();
			return result;
		}
		[Skip]
		public static Bernoulli SelectorAverageConditional(Bernoulli result)
		{
			result.SetToUniform();
			return result;
		}

		public static TDist ValueAverageConditional<TDist>([SkipIfUniform] IList<TDist> enterPartial, [SkipIfUniform] Discrete selector, TDist value, int[] indices, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>,
								SettableToRatio<TDist>, SettableToWeightedSum<TDist>, CanGetLogAverageOf<TDist>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (selector.Dimension < enterPartial.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				// TODO: use pre-allocated buffers
				double logProbSum = selector.GetLogProb(indices[0]);
				if (!double.IsNegativeInfinity(logProbSum)) {
					result.SetToProduct(value, enterPartial[0]);
				}
				if (indices.Length > 1) {
					TDist product = (TDist)value.Clone();
					for (int i = 1; i < indices.Length; i++) {
						double logProb = selector.GetLogProb(indices[i]);
						double shift = Math.Max(logProbSum, logProb);
						// avoid (-Infinity) - (-Infinity)
						if (Double.IsNegativeInfinity(shift)) {
							if (i == selector.Dimension - 1) {
								throw new AllZeroException();
							}
							// do nothing
						} else {
							double productWeight = Math.Exp(logProb - shift);
							if (productWeight > 0) {
								product.SetToProduct(value, enterPartial[i]);
								result.SetToSum(Math.Exp(logProbSum - shift), result, productWeight, product);
								logProbSum = MMath.LogSumExp(logProbSum, logProb);
							}
						}
					}
				}
				if (indices.Length < selector.Dimension) {
					double logProb = MMath.Log1MinusExp(logProbSum);
					double shift = Math.Max(logProbSum, logProb);
					if (Double.IsNegativeInfinity(shift)) throw new AllZeroException();
					result.SetToSum(Math.Exp(logProbSum - shift), result, Math.Exp(logProb - shift), value);
				}
				result.SetToRatio(result, value);
			}
			return result;
		}
		public static TDist ValueAverageConditional<TDist>(
			[SkipIfAllUniform] IList<TDist> enterPartial,
			int selector, int[] indices, TDist result)
			where TDist : IDistribution<T>, SettableTo<TDist>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				result.SetToUniform();
				for (int i = 0; i < indices.Length; i++) {
					if (selector == indices[i]) {
						result.SetTo(enterPartial[i]);
						break;
					}
				}
				return result;
			}
		}
#if false
		public static TDist ValueAverageConditional<TDist>(
			IList<T> enterPartial,
			int selector, int[] indices, TDist result)
			where TDist : IDistribution<T>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				result.SetToUniform();
				for (int i = 0; i < indices.Length; i++) {
					if (selector == indices[i]) {
						result.Point = enterPartial[i];
						break;
					}
				}
				return result;
			}
		}
#endif

		public static TDist ValueAverageConditional<TDist>([SkipIfUniform] IList<TDist> enterPartial, [SkipIfUniform] Bernoulli selector, TDist value, int[] indices, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>,
								SettableToRatio<TDist>, SettableToWeightedSum<TDist>, CanGetLogAverageOf<TDist>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (2 < enterPartial.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				// TODO: use pre-allocated buffers
				double logProbSum = (indices[0]==0) ? selector.GetLogProbTrue() : selector.GetLogProbFalse();
				if (!double.IsNegativeInfinity(logProbSum)) {
					result.SetToProduct(value, enterPartial[0]);
				}
				if (indices.Length > 1) {
					TDist product = (TDist)value.Clone();
					for (int i = 1; i < indices.Length; i++) {
						double logProb = (indices[i]==0) ? selector.GetLogProbTrue() : selector.GetLogProbFalse();
						double shift = Math.Max(logProbSum, logProb);
						// avoid (-Infinity) - (-Infinity)
						if (Double.IsNegativeInfinity(shift)) {
							if (i == 1) {
								throw new AllZeroException();
							}
							// do nothing
						} else {
							double productWeight = Math.Exp(logProb - shift);
							if (productWeight > 0) {
								product.SetToProduct(value, enterPartial[i]);
								result.SetToSum(Math.Exp(logProbSum - shift), result, productWeight, product);
								logProbSum = MMath.LogSumExp(logProbSum, logProb);
							}
						}
					}
				}
				if (indices.Length < 2) {
					double logProb = MMath.Log1MinusExp(logProbSum);
					double shift = Math.Max(logProbSum, logProb);
					if (Double.IsNegativeInfinity(shift)) throw new AllZeroException();
					result.SetToSum(Math.Exp(logProbSum - shift), result, Math.Exp(logProb - shift), value);
				}
				if (GateEnterOp<T>.ForceProper && (result is Gaussian)) {
					Gaussian r = (Gaussian)(object)result;
					r.SetToRatioProper(r, (Gaussian)(object)value);
					result = (TDist)(object)r;
				} else {
					result.SetToRatio(result, value);
				}
			}
			return result;
		}
		public static TDist ValueAverageConditional<TDist>(
			[SkipIfAllUniform] IList<TDist> enterPartial,
			bool selector, int[] indices, TDist result)
			where TDist : IDistribution<T>, SettableTo<TDist>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				int caseNumber = selector ? 0 : 1;
				result.SetToUniform();
				for (int i = 0; i < indices.Length; i++) {
					if (caseNumber == indices[i]) {
						result.SetTo(enterPartial[i]);
						break;
					}
				}
				return result;
			}
		}
		public static TDist ValueAverageConditional<TDist>(
			IList<T> enterPartial,
			bool selector, int[] indices, TDist result)
			where TDist : IDistribution<T>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				int caseNumber = selector ? 0 : 1;
				result.SetToUniform();
				for (int i = 0; i < indices.Length; i++) {
					if (caseNumber == indices[i]) {
						result.Point = enterPartial[i];
						break;
					}
				}
				return result;
			}
		}

#if false
		/// <summary>
		/// EP message to 'cases'
		/// </summary>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'cases' conditioned on the given values.
		/// </para></remarks>
		[Skip]
		public static BernoulliList CasesAverageConditional<BernoulliList>(BernoulliList result)
			where BernoulliList : SettableToUniform
		{
			result.SetToUniform();
			return result;
		}
		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="enterPartial">Incoming message from 'enterPartial'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="cases">Incoming message from 'cases'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(enterPartial,cases) p(enterPartial,cases) factor(enterPartial,cases,value,indices)]/p(value)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enterPartial"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="cases"/> is not a proper distribution</exception>
		public static TDist ValueAverageConditional<TDist>([SkipIfUniform] IList<TDist> enterPartial, [SkipIfUniform] IList<Bernoulli> cases, TDist value, int[] indices, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>,
								SettableToRatio<TDist>, SettableToWeightedSum<TDist>, CanGetLogAverageOf<TDist>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (cases.Count < enterPartial.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				// TODO: use pre-allocated buffers
				double logProbSum = cases[indices[0]].LogOdds;
				if (!double.IsNegativeInfinity(logProbSum)) {
					result.SetToProduct(value, enterPartial[0]);
				}
				if (indices.Length > 1) {
					TDist product = (TDist)value.Clone();
					for (int i = 1; i < indices.Length; i++) {
						double logProb = cases[indices[i]].LogOdds;
						double shift = Math.Max(logProbSum, logProb);
						// avoid (-Infinity) - (-Infinity)
						if (Double.IsNegativeInfinity(shift)) {
							if (i == cases.Count - 1) {
								throw new AllZeroException();
							}
							// do nothing
						} else {
							double productWeight = Math.Exp(logProb - shift);
							if (productWeight > 0) {
								product.SetToProduct(value, enterPartial[i]);
								result.SetToSum(Math.Exp(logProbSum - shift), result, productWeight, product);
								logProbSum = MMath.LogSumExp(logProbSum, logProb);
							}
						}
					}
				}
				if (indices.Length < cases.Count) {
					double logProb = MMath.Log1MinusExp(logProbSum);
					double shift = Math.Max(logProbSum, logProb);
					if (Double.IsNegativeInfinity(shift)) throw new AllZeroException();
					result.SetToSum(Math.Exp(logProbSum - shift), result, Math.Exp(logProb - shift), value);
				}
				if (GateEnterOp<T>.ForceProper && (result is Gaussian)) {
					Gaussian r = (Gaussian)(object)result;
					r.SetToRatioProper(r, (Gaussian)(object)value);
					result = (TDist)(object)r;
				} else {
					result.SetToRatio(result, value);
				}
			}
			return result;
		}
		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="enterPartial">Incoming message from 'enterPartial'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="cases">Constant value for 'cases'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(enterPartial,cases) p(enterPartial,cases) factor(enterPartial,cases,value,indices)]/p(value)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enterPartial"/> is not a proper distribution</exception>
		public static TDist ValueAverageConditional<TDist>(
			[SkipIfAllUniform] IList<TDist> enterPartial,
			IList<bool> cases, int[] indices, TDist result)
			where TDist : IDistribution<T>, SettableTo<TDist>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (cases.Count < enterPartial.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				result.SetToUniform();
				for (int i = 0; i < indices.Length; i++) {
					if (cases[indices[i]]) {
						result.SetTo(enterPartial[i]);
						break;
					}
				}
				return result;
			}
		}
		public static TDist ValueAverageConditional<TDist>(
			IList<T> enterPartial,
			IList<bool> cases, int[] indices, TDist result)
			where TDist : IDistribution<T>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (cases.Count < enterPartial.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				for (int i = 0; i < indices.Length; i++) {
					if (cases[indices[i]]) {
						result.Point = enterPartial[i];
						break;
					}
				}
				return result;
			}
		}
#endif

		//-- VMP ---------------------------------------------------------------------------------------------------------
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enterPartial,cases,value,indices))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }
		/// <summary>
		/// VMP message to 'enterPartial'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'enterPartial' as the random arguments are varied.
		/// The formula is <c>proj[sum_(value) p(value) factor(enterPartial,cases,value,indices)]</c>.
		/// </para></remarks>
		public static TList EnterPartialAverageLogarithm<T, TList>([IsReturnedInEveryElement] T value, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(value);
			return result;
		}

		public static TDist ValueAverageLogarithm<TDist>([SkipIfAllUniform] IList<TDist> enterPartial, [SkipIfUniform] Discrete selector, int[] indices, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToPower<TDist>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (selector.Dimension < enterPartial.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				double scale = selector[indices[0]];
				result.SetToPower(enterPartial[0], scale);
				if (indices.Length > 1) {
					// TODO: use pre-allocated buffer
					TDist power = (TDist)result.Clone();
					for (int i = 1; i < indices.Length; i++) {
						scale = selector[indices[i]];
						power.SetToPower(enterPartial[i], scale);
						result.SetToProduct(result, power);
					}
				}
			}
			return result;
		}
		public static TDist ValueAverageLogarithm<TDist>([SkipIfAllUniform] IList<TDist> enterPartial, [SkipIfUniform] Bernoulli selector, int[] indices, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToPower<TDist>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (2 < enterPartial.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				double scale = (indices[0]==0) ? selector.GetProbTrue() : selector.GetProbFalse();
				result.SetToPower(enterPartial[0], scale);
				if (indices.Length > 1) {
					// TODO: use pre-allocated buffer
					TDist power = (TDist)result.Clone();
					for (int i = 1; i < indices.Length; i++) {
						scale = (indices[i]==0) ? selector.GetProbTrue() : selector.GetProbFalse();
						power.SetToPower(enterPartial[i], scale);
						result.SetToProduct(result, power);
					}
				}
			}
			return result;
		}

		[Skip]
		public static Discrete SelectorAverageLogarithm(Discrete result)
		{
			result.SetToUniform();
			return result;
		}
		[Skip]
		public static Bernoulli SelectorAverageLogarithm(Bernoulli result)
		{
			result.SetToUniform();
			return result;
		}
		
#if false
		/// <summary>
		/// VMP message to 'cases'
		/// </summary>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'cases'.
		/// Because the factor is deterministic, 'enterPartial' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(value) p(value) log(sum_enterPartial p(enterPartial) factor(enterPartial,cases,value,indices)))</c>.
		/// </para></remarks>
		[Skip]
		public static BernoulliList CasesAverageLogarithm<BernoulliList>(BernoulliList result)
			where BernoulliList : SettableToUniform
		{
			result.SetToUniform();
			return result;
		}
		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="enterPartial">Incoming message from 'enterPartial'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="cases">Incoming message from 'cases'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'value'.
		/// Because the factor is deterministic, 'enterPartial' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(cases) p(cases) log(sum_enterPartial p(enterPartial) factor(enterPartial,cases,value,indices)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enterPartial"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="cases"/> is not a proper distribution</exception>
		public static TDist ValueAverageLogarithm<TDist>([SkipIfAllUniform] IList<TDist> enterPartial, [SkipIfUniform] IList<Bernoulli> cases, int[] indices, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToPower<TDist>
		{
			if (indices.Length != enterPartial.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (cases.Count < enterPartial.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				double scale = Math.Exp(cases[indices[0]].LogOdds);
				result.SetToPower(enterPartial[0], scale);
				if (indices.Length > 1) {
					// TODO: use pre-allocated buffer
					TDist power = (TDist)result.Clone();
					for (int i = 1; i < indices.Length; i++) {
						scale = Math.Exp(cases[indices[i]].LogOdds);
						power.SetToPower(enterPartial[i], scale);
						result.SetToProduct(result, power);
					}
				}
			}
			return result;
		}
#endif
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.EnterPartialTwo{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Gate), "EnterPartialTwo<>")]
	[Quality(QualityBand.Experimental)]
	public static class GateEnterPartialTwoOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enterPartialTwo,case0,case1,value,indices))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio()
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enterPartialTwo,case0,case1,value,indices))</c>.
		/// </para></remarks>
		[Skip]
		public static double LogAverageFactor()
		{
			return 0.0;
		}

		/// <summary>
		/// EP message to 'enterPartialTwo'
		/// </summary>
		/// <param name="value">Incoming message from 'value'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'enterPartialTwo' as the random arguments are varied.
		/// The formula is <c>proj[p(enterPartialTwo) sum_(value) p(value) factor(enterPartialTwo,case0,case1,value,indices)]/p(enterPartialTwo)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="value"/> is not a proper distribution</exception>
		public static TList EnterPartialTwoAverageConditional<T, TList>([SkipIfUniform] T value, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(value);
			return result;
		}
		/// <summary>
		/// EP message to 'case0'
		/// </summary>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'case0' conditioned on the given values.
		/// </para></remarks>
		[Skip]
		public static Bernoulli Case0AverageConditional(Bernoulli result)
		{
			result.SetToUniform();
			return result;
		}

		/// <summary>
		/// EP message to 'case1'
		/// </summary>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'case1' conditioned on the given values.
		/// </para></remarks>
		[Skip]
		public static Bernoulli Case1AverageConditional(Bernoulli result)
		{
			return Case0AverageConditional(result);
		}

		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="enterPartialTwo">Incoming message from 'enterPartialTwo'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="case0">Incoming message from 'case0'.</param>
		/// <param name="case1">Incoming message from 'case1'.</param>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(enterPartialTwo,case0,case1) p(enterPartialTwo,case0,case1) factor(enterPartialTwo,case0,case1,value,indices)]/p(value)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enterPartialTwo"/> is not a proper distribution</exception>
		public static T ValueAverageConditional<T>([SkipIfAllUniform] IList<T> enterPartialTwo, Bernoulli case0, Bernoulli case1, T value, int[] indices, T result)
			where T : ICloneable, SettableToUniform, SettableToProduct<T>,
								SettableToRatio<T>, SettableToWeightedSum<T>, CanGetLogAverageOf<T>
		{
			if (indices.Length != enterPartialTwo.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (2 < enterPartialTwo.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				// TODO: use pre-allocated buffers
				result.SetToProduct(value, enterPartialTwo[0]);
				double scale = Math.Exp((indices[0] ==0 ? case0 : case1).LogOdds);
				double sumCases = scale;
				double resultScale = scale;
				if (indices.Length > 1) {
					T product = (T)value.Clone();
					for (int i = 1; i < indices.Length; i++) {
						product.SetToProduct(value, enterPartialTwo[i]);
						scale = Math.Exp((indices[i] ==0 ? case0 : case1).LogOdds);
						result.SetToSum(resultScale, result, scale, product);
						resultScale += scale;
						sumCases += scale;
					}
				}
				double totalCases = Math.Exp(case0.LogOdds) + Math.Exp(case1.LogOdds);
				result.SetToSum(resultScale, result, totalCases - sumCases, value);
				result.SetToRatio(result, value);
			}
			return result;
		}

		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="enterPartialTwo">Incoming message from 'enterPartialTwo'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="case1">Constant value for 'case0'.</param>
		/// <param name="case2">Constant value for 'case1'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(enterPartialTwo,case0,case1) p(enterPartialTwo,case0,case1) factor(enterPartialTwo,case0,case1,value,indices)]/p(value)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enterPartialTwo"/> is not a proper distribution</exception>
		public static T ValueAverageConditional<T, TDomain>(
			[SkipIfAllUniform] IList<T> enterPartialTwo,
			bool case1, bool case2, int[] indices, T result)
			where T : IDistribution<TDomain>, SettableTo<T>
		{
			if (indices.Length != enterPartialTwo.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (2 < enterPartialTwo.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				result.SetToUniform();
				for (int i = 0; i < indices.Length; i++) {
					if ((indices[i] == 0 && case1) || (indices[i] == 1 && case2)) {
						result.SetTo(enterPartialTwo[indices[i]]);
						break;
					}
				}
				return result;
			}
		}

		//-- VMP ---------------------------------------------------------------------------------------------------------
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enterPartialTwo,case0,case1,value,indices))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }
		/// <summary>
		/// VMP message to 'enterPartialTwo'
		/// </summary>
		/// <param name="value">Incoming message from 'value'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'enterPartialTwo' as the random arguments are varied.
		/// The formula is <c>proj[sum_(value) p(value) factor(enterPartialTwo,case0,case1,value,indices)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="value"/> is not a proper distribution</exception>
		public static TList EnterPartialTwoAverageLogarithm<T, TList>([SkipIfUniform] T value, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(value);
			return result;
		}
		/// <summary>
		/// VMP message to 'case0'
		/// </summary>
		/// <param name="enterPartialTwo">Incoming message from 'enterPartialTwo'.</param>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'case0'.
		/// Because the factor is deterministic, 'enterPartialTwo' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(value) p(value) log(sum_enterPartialTwo p(enterPartialTwo) factor(enterPartialTwo,case0,case1,value,indices)))</c>.
		/// </para></remarks>
		[Skip]
		public static Bernoulli Case0AverageLogarithm<T>(IList<T> enterPartialTwo, T value, int[] indices, Bernoulli result)
		{
			result.SetToUniform();
			return result;
		}

		/// <summary>
		/// VMP message to 'case1'
		/// </summary>
		/// <param name="enterPartialTwo">Incoming message from 'enterPartialTwo'.</param>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'case1'.
		/// Because the factor is deterministic, 'enterPartialTwo' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(value) p(value) log(sum_enterPartialTwo p(enterPartialTwo) factor(enterPartialTwo,case0,case1,value,indices)))</c>.
		/// </para></remarks>
		[Skip]
		public static Bernoulli Case1AverageLogarithm<T>(IList<T> enterPartialTwo, T value, int[] indices, Bernoulli result)
		{
			result.SetToUniform();
			return result;
		}

		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="enterPartialTwo">Incoming message from 'enterPartialTwo'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="case0">Incoming message from 'case0'.</param>
		/// <param name="case1">Incoming message from 'case1'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'value'.
		/// Because the factor is deterministic, 'enterPartialTwo' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(case0,case1) p(case0,case1) log(sum_enterPartialTwo p(enterPartialTwo) factor(enterPartialTwo,case0,case1,value,indices)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enterPartialTwo"/> is not a proper distribution</exception>
		public static T ValueAverageLogarithm<T>([SkipIfAllUniform] IList<T> enterPartialTwo, Bernoulli case0, Bernoulli case1, int[] indices, T result)
	where T : ICloneable, SettableToProduct<T>, SettableToPower<T>
		{
			if (indices.Length != enterPartialTwo.Count) throw new ArgumentException("indices.Length != enterPartial.Count");
			if (2 < enterPartialTwo.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			if (indices.Length == 0) throw new ArgumentException("indices.Length == 0");
			else {
				double scale = Math.Exp((indices[0] ==0 ? case0 : case1).LogOdds);
				result.SetToPower(enterPartialTwo[0], scale);
				if (indices.Length > 1) {
					// TODO: use pre-allocated buffer
					T power = (T)result.Clone();
					for (int i = 1; i < indices.Length; i++) {
						scale = Math.Exp((indices[i] ==0 ? case0 : case1).LogOdds);
						power.SetToPower(enterPartialTwo[i], scale);
						result.SetToProduct(result, power);
					}
				}
			}
			return result;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.EnterOne{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Gate), "EnterOne<>", null, typeof(int), null, typeof(int))]
	[Quality(QualityBand.Mature)]
	public static class GateEnterOneOp<T>
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enterOne,cases,value,index))</c>.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio()
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enterOne,cases,value,index))</c>.
		/// </para></remarks>
		[Skip]
		public static double LogAverageFactor()
		{
			return 0.0;
		}

		/// <summary>
		/// EP message to 'enterOne'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'enterOne' as the random arguments are varied.
		/// The formula is <c>proj[p(enterOne) sum_(value) p(value) factor(enterOne,cases,value,index)]/p(enterOne)</c>.
		/// </para></remarks>
		public static T EnterOneAverageConditional<T>([IsReturned] T value)
		{
			return value;
		}

		[Skip]
		public static Discrete SelectorAverageConditional(Discrete result)
		{
			result.SetToUniform();
			return result;
		}
		public static TDist ValueAverageConditional<TDist>([SkipIfAllUniform] TDist enterOne, Discrete selector, TDist value, int index, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToRatio<TDist>, SettableToWeightedSum<TDist>, SettableTo<TDist>
		{
			double logProb = selector.GetLogProb(index);
			if (logProb == 0.0) {
				result.SetTo(enterOne);
			} else if (double.IsNegativeInfinity(logProb)) {
				result.SetToUniform();
			} else {
				result.SetToProduct(value, enterOne);
				double logOtherProb = MMath.Log1MinusExp(logProb);
				double shift = Math.Max(logProb, logOtherProb);
				// avoid (-Infinity) - (-Infinity)
				if (Double.IsNegativeInfinity(shift)) throw new AllZeroException();
				result.SetToSum(Math.Exp(logProb - shift), result, Math.Exp(logOtherProb - shift), value);
				if (GateEnterOp<T>.ForceProper && (result is Gaussian)) {
					Gaussian r = (Gaussian)(object)result;
					r.SetToRatioProper(r, (Gaussian)(object)value);
					result = (TDist)(object)r;
				} else {
					result.SetToRatio(result, value);
				}
			}
			return result;
		}
		public static TDist ValueAverageConditional<TDist>([SkipIfAllUniform] TDist enterOne, int selector, int index, TDist result)
			where TDist : IDistribution<T>, SettableTo<TDist>
		{
			if (selector == index)
				result.SetTo(enterOne);
			else
				result.SetToUniform();
			return result;
		}
	
#if false
		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		[Skip]
		public static Bernoulli BAverageConditional(Bernoulli result)
		{
			result.SetToUniform();
			return result;
		}
		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="enterOne">Incoming message from 'enterOne'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="value">Incoming message from 'value'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(enterOne,cases) p(enterOne,cases) factor(enterOne,cases,value,index)]/p(value)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enterOne"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="value"/> is not a proper distribution</exception>
		public static TDist ValueAverageConditional<TDist>([SkipIfAllUniform] TDist enterOne, Bernoulli b, [Proper] TDist value, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToRatio<TDist>, SettableToWeightedSum<TDist>, SettableTo<TDist>
		{
			double logProb = b.LogOdds;
			if (logProb == 0.0) {
				result.SetTo(enterOne);
			} else if (double.IsNegativeInfinity(logProb)) {
				result.SetToUniform();
			} else {
				result.SetToProduct(value, enterOne);
				double logOtherProb = MMath.Log1MinusExp(logProb);
				double shift = Math.Max(logProb, logOtherProb);
				// avoid (-Infinity) - (-Infinity)
				if (Double.IsNegativeInfinity(shift)) throw new AllZeroException();
				result.SetToSum(Math.Exp(logProb - shift), result, Math.Exp(logOtherProb - shift), value);
				if (GateEnterOp<T>.ForceProper && (result is Gaussian)) {
					Gaussian r = (Gaussian)(object)result;
					r.SetToRatioProper(r, (Gaussian)(object)value);
					result = (TDist)(object)r;
				} else {
					result.SetToRatio(result, value);
				}
			}
			return result;
		}
		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="enterOne">Incoming message from 'enterOne'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(enterOne,cases) p(enterOne,cases) factor(enterOne,cases,value,index)]/p(value)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enterOne"/> is not a proper distribution</exception>
		public static TDist ValueAverageConditional<TDist>([SkipIfAllUniform] TDist enterOne, bool b, TDist result)
			where TDist : IDistribution<T>, SettableTo<TDist>
		{
			if (b)
				result.SetTo(enterOne);
			else
				result.SetToUniform();
			return result;
		}
#endif

		//-- VMP ---------------------------------------------------------------------------------------------------------
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enterOne,cases,value,index))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'enterOne'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'enterOne' as the random arguments are varied.
		/// The formula is <c>proj[sum_(value) p(value) factor(enterOne,cases,value,index)]</c>.
		/// </para></remarks>
		public static T EnterOneAverageLogarithm<T>([IsReturned] T value)
		{
			return value;
		}

		[Skip]
		public static Discrete SelectorAverageLogarithm(Discrete result)
		{
			result.SetToUniform();
			return result;
		}
		public static TDist ValueAverageLogarithm<TDist>([SkipIfUniform] TDist enterOne, Discrete selector, int index, TDist result)
			where TDist : SettableToPower<TDist>
		{
			double scale = selector[index];
			result.SetToPower(enterOne, scale);
			return result;
		}
	
#if false
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		[Skip]
		public static Bernoulli BAverageLogarithm(Bernoulli result)
		{
			result.SetToUniform();
			return result;
		}
		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="enterOne">Incoming message from 'enterOne'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		public static TDist ValueAverageLogarithm<TDist>([SkipIfUniform] TDist enterOne, Bernoulli b, TDist result)
			where TDist : SettableToPower<TDist>
		{
			double scale = Math.Exp(b.LogOdds);
			result.SetToPower(enterOne, scale);
			return result;
		}
#endif
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.Enter{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Gate), "Enter<>", null, typeof(bool), null)]
	[FactorMethod(typeof(Gate), "Enter<>", null, typeof(int), null)]
	[Quality(QualityBand.Stable)]
	public static class GateEnterOp<T>
	{
		/// <summary>
		/// Force proper messages
		/// </summary>
		public static bool ForceProper;

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enter,cases,value))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio()
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enter,cases,value))</c>.
		/// </para></remarks>
		[Skip]
		public static double LogAverageFactor()
		{
			return 0.0;
		}

		/// <summary>
		/// EP message to 'enter'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'enter' as the random arguments are varied.
		/// The formula is <c>proj[p(enter) sum_(value) p(value) factor(enter,cases,value)]/p(enter)</c>.
		/// </para></remarks>
		public static TList EnterAverageConditional<T, TList>([IsReturnedInEveryElement] T value, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(value);
			return result;
		}
		[Skip]
		public static ArrayType EnterInit<T, ArrayType>(Discrete selector, [IgnoreDependency] T value, IArrayFactory<T, ArrayType> factory)
			where T : ICloneable
		{
			return factory.CreateArray(selector.Dimension, i => (T)value.Clone());
		}

		[Skip]
		public static Discrete SelectorAverageConditional(Discrete result)
		{
			result.SetToUniform();
			return result;
		}

		public static TDist ValueAverageConditional<TDist>([SkipIfAllUniform] IList<TDist> enter, Discrete selector, TDist value, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>,
			SettableToRatio<TDist>, SettableToWeightedSum<TDist>, CanGetLogAverageOf<TDist>
		{
			if (selector.Dimension != enter.Count) throw new ArgumentException("selector.Dimension != enter.Count");
			// TODO: use pre-allocated buffers
			double logProbSum = selector.GetLogProb(0);
			if (!double.IsNegativeInfinity(logProbSum)) {
				result.SetToProduct(value, enter[0]);
			}
			if (selector.Dimension > 1) {
				TDist product = (TDist)value.Clone();
				for (int i = 1; i < selector.Dimension; i++) {
					double logProb = selector.GetLogProb(i);
					double shift = Math.Max(logProbSum, logProb);
					// avoid (-Infinity) - (-Infinity)
					if (Double.IsNegativeInfinity(shift)) {
						if (i == selector.Dimension - 1) {
							throw new AllZeroException();
						}
						// do nothing
					} else {
						double productWeight = Math.Exp(logProb - shift);
						if (productWeight > 0) {
							product.SetToProduct(value, enter[i]);
							result.SetToSum(Math.Exp(logProbSum - shift), result, productWeight, product);
							logProbSum = MMath.LogSumExp(logProbSum, logProb);
						}
					}
				}
			}
			result.SetToRatio(result, value);
			return result;
		}
		public static TDist ValueAverageConditional<TDist>([SkipIfAllUniform] IList<TDist> enter, int selector, TDist result)
			where TDist : IDistribution<T>, SettableTo<TDist>
		{
			result.SetTo(enter[selector]);
			return result;
		}

#if false
		/// <summary>
		/// EP message to 'cases'
		/// </summary>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'cases' conditioned on the given values.
		/// </para></remarks>
		[Skip]
		public static BernoulliList CasesAverageConditional<BernoulliList>(BernoulliList result)
			where BernoulliList : SettableToUniform
		{
			result.SetToUniform();
			return result;
		}
		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="enter">Incoming message from 'enter'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(enter,cases) p(enter,cases) factor(enter,cases,value)]/p(value)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enter"/> is not a proper distribution</exception>
		public static TDist ValueAverageConditional<TDist>([SkipIfAllUniform] IList<TDist> enter, IList<Bernoulli> cases, TDist value, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>,
			SettableToRatio<TDist>, SettableToWeightedSum<TDist>, CanGetLogAverageOf<TDist>
		{
			if (cases.Count < enter.Count) throw new ArgumentException("cases.Count < enter.Count");
			// TODO: use pre-allocated buffers
			double logProbSum = cases[0].LogOdds;
			if (!double.IsNegativeInfinity(logProbSum)) {
				result.SetToProduct(value, enter[0]);
			}
			if (cases.Count > 1) {
				TDist product = (TDist)value.Clone();
				for (int i = 1; i < cases.Count; i++) {
					double logProb = cases[i].LogOdds;
					double shift = Math.Max(logProbSum, logProb);
					// avoid (-Infinity) - (-Infinity)
					if (Double.IsNegativeInfinity(shift)) {
						if (i == cases.Count - 1) {
							throw new AllZeroException();
						}
						// do nothing
					} else {
						double productWeight = Math.Exp(logProb - shift);
						if (productWeight > 0) {
							product.SetToProduct(value, enter[i]);
							result.SetToSum(Math.Exp(logProbSum - shift), result, productWeight, product);
							logProbSum = MMath.LogSumExp(logProbSum, logProb);
						}
					}
				}
			}
			if (ForceProper && (result is Gaussian)) {
				Gaussian r = (Gaussian)(object)result;
				r.SetToRatioProper(r, (Gaussian)(object)value);
				result = (TDist)(object)r;
			} else {
				result.SetToRatio(result, value);
			}
			return result;
		}
		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="enter">Incoming message from 'enter'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="cases">Constant value for 'cases'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(enter) p(enter) factor(enter,cases,value)]/p(value)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enter"/> is not a proper distribution</exception>
		public static TDist ValueAverageConditional<TDist>([SkipIfAllUniform] IList<TDist> enter, bool[] cases, TDist result)
			where TDist : IDistribution<T>, SettableTo<TDist>
		{
			if (cases.Length < enter.Count) throw new ArgumentException("cases.Count < enter.Count");
			result.SetToUniform();
			for (int i = 0; i < cases.Length; i++) {
				if (cases[i]) {
					result.SetTo(enter[i]);
					break;
				}
			}
			return result;
		}
#endif

		//-- VMP ---------------------------------------------------------------------------------------------------------
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(enter,cases,value))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'enter'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'enter' as the random arguments are varied.
		/// The formula is <c>proj[sum_(value) p(value) factor(enter,cases,value)]</c>.
		/// </para></remarks>
		public static TList EnterAverageLogarithm<T, TList>([IsReturnedInEveryElement] T value, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(value);
			return result;
		}

		[Skip]
		public static Discrete SelectorAverageLogarithm(Discrete result)
		{
			result.SetToUniform();
			return result;
		}
		public static TDist ValueAverageLogarithm<TDist>([SkipIfAllUniform] IList<TDist> enter, Discrete selector, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToPower<TDist>
		{
			if (selector.Dimension != enter.Count) throw new ArgumentException("selector.Dimension != enterPartial.Count");
			double scale = selector[0];
			result.SetToPower(enter[0], scale);
			if (selector.Dimension > 1) {
				// TODO: use pre-allocated buffer
				TDist power = (TDist)result.Clone();
				for (int i = 1; i < selector.Dimension; i++) {
					scale = selector[i];
					power.SetToPower(enter[i], scale);
					result.SetToProduct(result, power);
				}
			}
			return result;
		}
	
#if false
		/// <summary>
		/// VMP message to 'cases'
		/// </summary>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'cases' conditioned on the given values.
		/// </para></remarks>
		[Skip]
		public static BernoulliList CasesAverageLogarithm<BernoulliList>(BernoulliList result)
			where BernoulliList : SettableToUniform
		{
			result.SetToUniform();
			return result;
		}
		// result = prod_i enterPartial[i]^cases[indices[i]]
		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="enter">Incoming message from 'enter'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'value'.
		/// Because the factor is deterministic, 'enter' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(cases) p(cases) log(sum_enter p(enter) factor(enter,cases,value)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="enter"/> is not a proper distribution</exception>
		public static TDist ValueAverageLogarithm<TDist>([SkipIfAllUniform] IList<TDist> enter, IList<Bernoulli> cases, TDist result)
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToPower<TDist>
		{
			if (cases.Count < enter.Count) throw new ArgumentException("cases.Count < enterPartial.Count");
			double scale = Math.Exp(cases[0].LogOdds);
			result.SetToPower(enter[0], scale);
			if (cases.Count > 1) {
				// TODO: use pre-allocated buffer
				TDist power = (TDist)result.Clone();
				for (int i = 1; i < cases.Count; i++) {
					scale = Math.Exp(cases[i].LogOdds);
					power.SetToPower(enter[i], scale);
					result.SetToProduct(result, power);
				}
			}
			return result;
		}
#endif
	}

}
