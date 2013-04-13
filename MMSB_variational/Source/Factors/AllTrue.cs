// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.AllTrue"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "AllTrue")]
	[Quality(QualityBand.Preview)]
	public static class AllTrueOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="allTrue">Constant value for 'allTrue'.</param>
		/// <param name="array">Constant value for 'array'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(allTrue,array))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool allTrue, IList<bool> array)
		{
			return (allTrue == Factor.AllTrue(array)) ? 0.0 : Double.NegativeInfinity;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="allTrue">Constant value for 'allTrue'.</param>
		/// <param name="array">Constant value for 'array'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(allTrue,array))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool allTrue, IList<bool> array) { return LogAverageFactor(allTrue, array); }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="allTrue">Constant value for 'allTrue'.</param>
		/// <param name="array">Constant value for 'array'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(allTrue,array))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(bool allTrue, IList<bool> array) { return LogAverageFactor(allTrue, array); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="allTrue">Incoming message from 'allTrue'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_allTrue">Outgoing message to 'allTrue'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(allTrue) p(allTrue) factor(allTrue,array))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="allTrue"/> is not a proper distribution</exception>
		public static double LogAverageFactor([SkipIfUniform] Bernoulli allTrue, [Fresh] Bernoulli to_allTrue)
		{
			return to_allTrue.GetLogAverageOf(allTrue);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="allTrue">Constant value for 'allTrue'.</param>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array) p(array) factor(allTrue,array))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static double LogAverageFactor(bool allTrue, [SkipIfAnyUniform] IList<Bernoulli> array)
		{
			Bernoulli to_allTrue = AllTrueAverageConditional(array);
			return to_allTrue.GetLogProb(allTrue);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="allTrue">Incoming message from 'allTrue'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks>returns 0</remarks>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli allTrue) { return 0.0; }

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="allTrue">Constant value for 'allTrue'.</param>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array) p(array) factor(allTrue,array))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio(bool allTrue, [SkipIfAnyUniform] IList<Bernoulli> array)
		{
			return LogAverageFactor(allTrue, array);
		}

		/// <summary>
		/// EP message to 'allTrue'.
		/// </summary>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <returns>The outgoing EP message to the 'allTrue' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'allTrue' as the random arguments are varied.
		/// The formula is <c>proj[p(allTrue) sum_(array) p(array) factor(allTrue,array)]/p(allTrue)</c>.
		/// </para></remarks>
		public static Bernoulli AllTrueAverageConditional(IList<Bernoulli> array)
		{
			double logOdds = Double.NegativeInfinity;
			for (int i = 0; i < array.Count; i++) {
				logOdds = Bernoulli.Or(logOdds, -array[i].LogOdds);
			}
			return Bernoulli.FromLogOdds(-logOdds);
		}

		/// <summary>
		/// EP message to 'allTrue'.
		/// </summary>
		/// <param name="array">Constant value for 'array'.</param>
		/// <returns>The outgoing EP message to the 'allTrue' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'allTrue' conditioned on the given values.
		/// </para></remarks>
		public static Bernoulli AllTrueAverageConditional(IList<bool> array)
		{
			foreach (bool b in array)
				if (!b)
					return Bernoulli.PointMass(false);
			return Bernoulli.PointMass(true);
		}

		/// <summary>
		/// EP message to 'array'.
		/// </summary>
		/// <param name="allTrue">Incoming message from 'allTrue'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'array' as the random arguments are varied.
		/// The formula is <c>proj[p(array) sum_(allTrue) p(allTrue) factor(allTrue,array)]/p(array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="allTrue"/> is not a proper distribution</exception>
		public static BernoulliList ArrayAverageConditional<BernoulliList>([SkipIfUniform] Bernoulli allTrue, IList<Bernoulli> array, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
		{
			if (result.Count == 0) { } else if (result.Count == 1) result[0] = allTrue;
			else if (result.Count == 2) {
				result[0] = BooleanAndOp.AAverageConditional(allTrue, array[1]);
				result[1] = BooleanAndOp.BAverageConditional(allTrue, array[0]);
			} else { // result.Count >= 3
				double notallTruePrevious = Double.NegativeInfinity;
				double[] notallTrueNext = new double[result.Count];
				notallTrueNext[notallTrueNext.Length-1] = Double.NegativeInfinity;
				for (int i = notallTrueNext.Length - 2; i >= 0; i--) {
					notallTrueNext[i] = Bernoulli.Or(-array[i+1].LogOdds, notallTrueNext[i+1]);
				}
				for (int i = 0; i < result.Count; i++) {
					double notallTrueExcept = Bernoulli.Or(notallTruePrevious, notallTrueNext[i]);
					result[i] = Bernoulli.FromLogOdds(-Bernoulli.Gate(-allTrue.LogOdds, notallTrueExcept));
					notallTruePrevious = Bernoulli.Or(notallTruePrevious, -array[i].LogOdds);
				}
			}
			return result;
		}

		/// <summary>
		/// EP message to 'array'.
		/// </summary>
		/// <param name="allTrue">Constant value for 'allTrue'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' conditioned on the given values.
		/// </para></remarks>
		public static BernoulliList ArrayAverageConditional<BernoulliList>(bool allTrue, IList<Bernoulli> array, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
		{
			return ArrayAverageConditional(Bernoulli.PointMass(allTrue), array, result);
		}

		/// <summary>
		/// EP message to 'array'.
		/// </summary>
		/// <param name="allTrue">Incoming message from 'allTrue'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'array' as the random arguments are varied.
		/// The formula is <c>proj[p(array) sum_(allTrue) p(allTrue) factor(allTrue,array)]/p(array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="allTrue"/> is not a proper distribution</exception>
		public static BernoulliList ArrayAverageConditional<BernoulliList>([SkipIfUniform] Bernoulli allTrue, IList<bool> array, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
		{
			if (result.Count == 0) { } else if (result.Count == 1) result[0] = allTrue;
			else if (result.Count == 2) {
				result[0] = BooleanAndOp.AAverageConditional(allTrue, array[1]);
				result[1] = BooleanAndOp.BAverageConditional(allTrue, array[0]);
			} else { // result.Count >= 3
				int trueCount = 0;
				int firstFalseIndex = -1;
				for (int i = 0; i < array.Count; i++) {
					if (array[i])
						trueCount++;
					else if (firstFalseIndex < 0)
						firstFalseIndex = i;
				}
				if (trueCount == array.Count) {
					for (int i=0; i < result.Count; i++)
						result[i] = BooleanAndOp.AAverageConditional(allTrue, true);
				} else {
					for (int i=0; i < result.Count; i++)
						result[i] = BooleanAndOp.AAverageConditional(allTrue, false);
					if (trueCount == result.Count - 1)
						result[firstFalseIndex] = BooleanAndOp.AAverageConditional(allTrue, true);
				}
			}
			return result;
		}

		// VMP //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		const string NotSupportedMessage = "Variational Message Passing does not support an AllTrue factor with fixed output.";

		/// <summary>
		/// VMP message to 'allTrue'.
		/// </summary>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <returns>The outgoing VMP message to the 'allTrue' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'allTrue' as the random arguments are varied.
		/// The formula is <c>proj[sum_(array) p(array) factor(allTrue,array)]</c>.
		/// </para></remarks>
		public static Bernoulli AllTrueAverageLogarithm(IList<Bernoulli> array)
		{
			return AllTrueAverageConditional(array);
		}

		/// <summary>
		/// VMP message to 'array'.
		/// </summary>
		/// <param name="allTrue">Incoming message from 'allTrue'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' with 'allTrue' integrated out.
		/// The formula is <c>sum_allTrue p(allTrue) factor(allTrue,array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="allTrue"/> is not a proper distribution</exception>
		public static BernoulliList ArrayAverageLogarithm<BernoulliList>([SkipIfUniform] Bernoulli allTrue, IList<Bernoulli> array, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
		{
			// when 'allTrue' is marginalized, the factor is proportional to exp(allTrue.LogOdds*prod_i a[i])
			// therefore we maintain the value of prod_i E(a[i]) as we update each a[i]'s distribution
			double prodProbTrue = 1.0;
			for (int i = 0; i < array.Count; i++) {
				prodProbTrue *= array[i].GetProbTrue();
			}
			for (int i = 0; i < result.Count; i++) {
				prodProbTrue /= array[i].GetProbTrue();
				Bernoulli ratio = array[i] / result[i];
				result[i] = Bernoulli.FromLogOdds(allTrue.LogOdds * prodProbTrue);
				Bernoulli newMarginal = ratio * result[i];
				prodProbTrue *= newMarginal.GetProbTrue();
			}
			return result;
		}

		/// <summary>
		/// VMP message to 'array'.
		/// </summary>
		/// <param name="allTrue">Constant value for 'allTrue'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' conditioned on the given values.
		/// </para></remarks>
		[NotSupported(AllTrueOp.NotSupportedMessage)]
		public static BernoulliList ArrayAverageLogarithm<BernoulliList>(bool allTrue, [Stochastic] IList<Bernoulli> array, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
		{
			throw new NotSupportedException(NotSupportedMessage);
			//return ArrayAverageLogarithm(Bernoulli.PointMass(allTrue), array, result);
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }
	}
}
