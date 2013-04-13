// (C) Copyright 2008 Microsoft Research Cambridge
#define UseRatioDir
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
#if true
	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.ExitingVariable{T}"/>, given random arguments to the function.
	/// </summary>
	/// <remarks><para>
	/// This factor is like ReplicateWithMarginal except Uses[0] plays the role of Def, and Def is
	/// considered a Use.  Needed only when a variable exits a gate in VMP.
	/// </para></remarks>
	[FactorMethod(typeof(Gate), "ExitingVariable<>")]
	[Quality(QualityBand.Preview)]
	public static class ExitingVariableOp
	{
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		public static T MarginalAverageLogarithm<T>([IsReturned] T Use)
		{
			return Use;
		}
		[Skip]
		public static T MarginalAverageLogarithmInit<T>(T Def)
			where T : ICloneable
		{
			return (T)Def.Clone();
		}

		public static T UseAverageLogarithm<T>([IsReturned] T Def)
		{
			return Def;
		}

		public static T DefAverageLogarithm<T>([IsReturned] T Use)
		{
			return Use;
		}
	}
	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.ReplicateExiting{T}"/>, given random arguments to the function.
	/// </summary>
	/// <remarks><para>
	/// This factor is like Replicate except Uses[0] plays the role of Def, and Def is
	/// considered a Use.  Needed only when a variable exits a gate in VMP.
	/// </para></remarks>
	[FactorMethod(typeof(Gate), "ReplicateExiting<>")]
	[Quality(QualityBand.Preview)]
	public static class ReplicateExitingOp
	{
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'Uses'.
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'Uses'.
		/// The formula is <c>int log(f(Uses,x)) q(x) dx</c> where <c>x = (Def,Marginal)</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		public static T UsesAverageLogarithm<T>([AllExceptIndex] IList<T> Uses, T Def, int resultIndex, T result)
				where T : SettableTo<T>, SettableToProduct<T>
		{
			if (resultIndex == 0) {
				result.SetTo(Def);
				result = Distribution.SetToProductWithAllExcept(result, Uses, 0);
			} else {
				result.SetTo(Uses[0]);
			}
			return result;
		}
		[Skip]
		public static T UsesAverageLogarithmInit<T>(T Def, int resultIndex)
			where T : ICloneable
		{
			return (T)Def.Clone();
		}

		/// <summary>
		/// VMP message to 'Def'.
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'Def'.
		/// The formula is <c>int log(f(Def,x)) q(x) dx</c> where <c>x = (Uses,Marginal)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		public static T DefAverageLogarithm<T>([SkipIfAllUniform] IList<T> Uses, T result)
		where T : SettableTo<T>
		{
			result.SetTo(Uses[0]);
			return result;
		}
	}
#else
	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.ExitingVariable{T}"/>, given random arguments to the function.
	/// </summary>
	/// <remarks><para>
	/// This factor is like ReplicateWithMarginal except Uses[0] plays the role of Def, and Def is
	/// considered a Use.  Needed only when a variable exits a gate in VMP.
	/// </para></remarks>
	[FactorMethod(typeof(Gate), "ExitingVariable<>")]
	[Quality(QualityBand.Preview)]
	public static class ExitingVariableOp
	{
		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (Uses,Def,Marginal)</c>.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'Marginal'.
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'Marginal'.
		/// The formula is <c>int log(f(Marginal,x)) q(x) dx</c> where <c>x = (Uses,Def)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		public static T MarginalAverageLogarithm<T>([SkipIfAllUniform] IList<T> Uses, T result)
		where T : SettableTo<T>
		{
			result.SetTo(Uses[0]);
			return result;
		}

		/// <summary>
		/// VMP message to 'Uses'.
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'Uses'.
		/// The formula is <c>int log(f(Uses,x)) q(x) dx</c> where <c>x = (Def,Marginal)</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		public static T UsesAverageLogarithm<T>([AllExceptIndex] IList<T> Uses, T Def, int resultIndex, T result)
				where T : SettableTo<T>, SettableToProduct<T>
		{
			if (resultIndex == 0) {
				result.SetTo(Def);
				result = Distribution.SetToProductWithAllExcept(result, Uses, 0);
			} else {
				result.SetTo(Uses[0]);
			}
			return result;
		}

		/// <summary>
		/// VMP message to 'Def'.
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'Def'.
		/// The formula is <c>int log(f(Def,x)) q(x) dx</c> where <c>x = (Uses,Marginal)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		public static T DefAverageLogarithm<T>([SkipIfAllUniform] IList<T> Uses, T result)
		where T : SettableTo<T>
		{
			return MarginalAverageLogarithm(Uses, result);
		}
	}
#endif

	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.Exit{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Gate), "Exit<>")]
	[Quality(QualityBand.Mature)]
	public static class GateExitOp<T>
	{
		public static bool ForceProper;

		[Skip]
		public static double LogEvidenceRatio<TDist>(TDist exit, IList<bool> cases, IList<T> values)
			where TDist : IDistribution<T>, CanGetLogAverageOf<TDist> 
		{ 
			return 0.0; 
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="exit">Incoming message from 'exit'.</param>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="values">Incoming message from 'values'.</param>
		/// <param name="to_exit">Outgoing message to 'exit'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(exit,cases,values) p(exit,cases,values) factor(exit,cases,values) / sum_exit p(exit) messageTo(exit))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<TDist>([SkipIfUniform] TDist exit, IList<Bernoulli> cases, IList<TDist> values, [Fresh] TDist to_exit)
			where TDist : IDistribution<T>, CanGetLogAverageOf<TDist>
		{
			return -to_exit.GetLogAverageOf(exit);
		}

		/// <summary>
		/// EP message to 'values'
		/// </summary>
		/// <param name="exit">Incoming message from 'exit'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'values' as the random arguments are varied.
		/// The formula is <c>proj[p(values) sum_(exit) p(exit) factor(exit,cases,values)]/p(values)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="exit"/> is not a proper distribution</exception>
		public static TList ValuesAverageConditional<T, TList>([IsReturnedInEveryElement] T exit, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(exit);
			return result;
		}
		/// <summary>
		/// EP message to 'cases'
		/// </summary>
		/// <param name="exit">Incoming message from 'exit'.</param>
		/// <param name="values">Incoming message from 'values'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'cases' as the random arguments are varied.
		/// The formula is <c>proj[p(cases) sum_(exit,values) p(exit,values) factor(exit,cases,values)]/p(cases)</c>.
		/// </para></remarks>
		public static BernoulliList CasesAverageConditional<TDist, BernoulliList>([SkipIfUniform] TDist exit, IList<TDist> values, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
			where TDist : IDistribution<T>, CanGetLogAverageOf<TDist>
		{
			for (int i = 0; i < values.Count; i++) {
				result[i] = Bernoulli.FromLogOdds(exit.GetLogAverageOf(values[i]));
			}
			return result;
		}
		/// <summary>
		/// EP message to 'exit'
		/// </summary>
		/// <param name="exit">Incoming message from 'exit'.</param>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'exit' as the random arguments are varied.
		/// The formula is <c>proj[p(exit) sum_(cases,values) p(cases,values) factor(exit,cases,values)]/p(exit)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		public static TDist ExitAverageConditional<TDist>(TDist exit, IList<Bernoulli> cases, [SkipIfUniform] IList<TDist> values, TDist result)
			where TDist : IDistribution<T>, SettableTo<TDist>, SettableToProduct<TDist>,
			SettableToRatio<TDist>, SettableToWeightedSum<TDist>, CanGetLogAverageOf<TDist>
		{
			if (cases.Count != values.Count) throw new ArgumentException("cases.Count != values.Count");
			if (cases.Count == 0) throw new ArgumentException("cases.Count == 0");
			else if (cases.Count == 1) {
				result.SetTo(values[0]);
			} else {
				double resultScale = Math.Exp(exit.GetLogAverageOf(values[0]) + cases[0].LogOdds);
				if (double.IsNaN(resultScale)) throw new AllZeroException();
				if (resultScale > 0) {
					result.SetToProduct(exit, values[0]);
				}
				// TODO: use pre-allocated buffer
				TDist product = (TDist)exit.Clone();
				for (int i = 1; i < cases.Count; i++) {
					double scale = Math.Exp(exit.GetLogAverageOf(values[i]) + cases[i].LogOdds);
					if (scale > 0) {
						product.SetToProduct(exit, values[i]);
						result.SetToSum(resultScale, result, scale, product);
						resultScale += scale;
					}
				}
				if (ForceProper && (result is Gaussian)) {
					Gaussian r = (Gaussian)(object)result;
					r.SetToRatioProper(r, (Gaussian)(object)exit);
					result = (TDist)(object)r;
				} else {
					result.SetToRatio(result, exit);
				}
			}
			return result;
		}
		/// <summary>
		/// EP message to 'exit'.
		/// </summary>
		/// <param name="exit">Incoming message from 'exit'.</param>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'exit'.
		/// The formula is <c>int f(exit,x) q(x) dx</c> where <c>x = (cases,values)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		public static TDist ExitAverageConditional2<TDist>(TDist exit, IList<Bernoulli> cases, [SkipIfAllUniform] IList<TDist> values, TDist result)
			where TDist : SettableTo<TDist>, ICloneable, SettableToProduct<TDist>,
			SettableToRatio<TDist>, SettableToWeightedSum<TDist>
		{
			if (cases.Count != values.Count) throw new ArgumentException("cases.Count != values.Count");
			if (cases.Count == 0) throw new ArgumentException("cases.Count == 0");
			else if (cases.Count == 1) {
				result.SetTo(values[0]);
			} else {
				result.SetToProduct(exit, values[0]);
				double scale = cases[0].LogOdds;
				double resultScale = scale;
				// TODO: use pre-allocated buffer
				TDist product = (TDist)exit.Clone();
				for (int i = 1; i < cases.Count; i++) {
					product.SetToProduct(exit, values[i]);
					scale = cases[i].LogOdds;
					double shift = Math.Max(resultScale, scale);
					// avoid (-Infinity) - (-Infinity)
					if (Double.IsNegativeInfinity(shift)) {
						if (i == cases.Count - 1) {
							throw new AllZeroException();
						}
						// do nothing
					} else {
						result.SetToSum(Math.Exp(resultScale - shift), result, Math.Exp(scale - shift), product);
						resultScale = MMath.LogSumExp(resultScale, scale);
					}
				}
				if (ForceProper && (result is Gaussian)) {
					Gaussian r = (Gaussian)(object)result;
					r.SetToRatioProper(r, (Gaussian)(object)exit);
					result = (TDist)(object)r;
				} else {
					result.SetToRatio(result, exit);
				}
			}
			return result;
		}
		/// <summary>
		/// EP message to 'cases'
		/// </summary>
		/// <param name="exit">Incoming message from 'exit'.</param>
		/// <param name="values">Incoming message from 'values'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'cases' as the random arguments are varied.
		/// The formula is <c>proj[p(cases) sum_(exit,values) p(exit,values) factor(exit,cases,values)]/p(cases)</c>.
		/// </para></remarks>
		public static BernoulliList CasesAverageConditional<TDist, BernoulliList>(TDist exit, IList<T> values, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
			where TDist : CanGetLogProb<T>
		{
			for (int i = 0; i < values.Count; i++) {
				result[i] = Bernoulli.FromLogOdds(exit.GetLogProb(values[i]));
			}
			return result;
		}
		/// <summary>
		/// EP message to 'exit'
		/// </summary>
		/// <param name="exit">Incoming message from 'exit'.</param>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="values">Incoming message from 'values'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'exit' as the random arguments are varied.
		/// The formula is <c>proj[p(exit) sum_(cases,values) p(cases,values) factor(exit,cases,values)]/p(exit)</c>.
		/// </para></remarks>
		public static TDist ExitAverageConditional<TDist>(
			TDist exit, IList<Bernoulli> cases, IList<T> values, TDist result)
			where TDist : ICloneable, HasPoint<T>, SettableToWeightedSum<TDist>
		{
			if (cases.Count != values.Count) throw new ArgumentException("cases.Count != values.Count");
			if (cases.Count == 0) throw new ArgumentException("cases.Count == 0");
			else if (cases.Count == 1) {
				result.Point = values[0];
			} else {
				result.Point = values[0];
				double scale = cases[0].LogOdds;
				double resultScale = scale;
				// TODO: overload SetToSum to accept constants.
				TDist product = (TDist)exit.Clone();
				for (int i = 1; i < cases.Count; i++) {
					product.Point = values[i];
					scale = cases[i].LogOdds;
					double shift = Math.Max(resultScale, scale);
					// avoid (-Infinity) - (-Infinity)
					if (Double.IsNegativeInfinity(shift)) {
						if (i == cases.Count - 1) {
							throw new AllZeroException();
						}
						// do nothing
					} else {
						result.SetToSum(Math.Exp(resultScale - shift), result, Math.Exp(scale - shift), product);
						resultScale = MMath.LogSumExp(resultScale, scale);
					}
				}
			}
			return result;
		}

		[Skip]
		public static TDist ExitAverageConditionalInit<TDist>([IgnoreDependency] IList<TDist> values)
			where TDist : ICloneable
		{
			return (TDist)values[0].Clone();
		}

		//-- VMP ------------------------------------------------------------------------------------------------
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="exit">Incoming message from 'exit'.</param>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_exit">Outgoing message to 'exit'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		public static double AverageLogFactor<TDist>([SkipIfUniform] TDist exit, [Fresh] TDist to_exit)
			where TDist : IDistribution<T>, CanGetAverageLog<TDist>
		{
			// cancel the evidence message from the child variable's child factors
			return -to_exit.GetAverageLog(exit);
		}
		/// <summary>
		/// VMP message to 'values'
		/// </summary>
		/// <param name="exit">Incoming message from 'exit'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'values' with 'exit' integrated out.
		/// The formula is <c>sum_exit p(exit) factor(exit,cases,values)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="exit"/> is not a proper distribution</exception>
		public static TList ValuesAverageLogarithm<T, TList>([IsReturnedInEveryElement] T exit, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(exit);
			return result;
		}
#if true
		/// <summary>
		/// VMP message to 'cases'
		/// </summary>
		/// <param name="exit">Incoming message from 'exit'.</param>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'cases'.
		/// Because the factor is deterministic, 'exit' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(values) p(values) log(sum_exit p(exit) factor(exit,cases,values)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		[NoTriggers] // see VmpTests.GateExitTriggerTest
		public static BernoulliList CasesAverageLogarithm<TDist, BernoulliList>([SkipIfUniform] TDist exit, [SkipIfAllUniform, Proper, Trigger] IList<TDist> values, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
			where TDist : CanGetAverageLog<TDist>
		{
			for (int i = 0; i < values.Count; i++) {
				result[i] = Bernoulli.FromLogOdds(values[i].GetAverageLog(exit));
			}
			return result;
		}
		// result = prod_i values[i]^cases[i]  (messages out of a gate are blurred)
		/// <summary>
		/// VMP message to 'exit'
		/// </summary>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'exit' as the random arguments are varied.
		/// The formula is <c>proj[sum_(cases,values) p(cases,values) factor(exit,cases,values)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		public static TDist ExitAverageLogarithm<TDist>(IList<Bernoulli> cases, [SkipIfAllUniform, Proper] IList<TDist> values, TDist result)
			where TDist : ICloneable, SettableToProduct<TDist>,
								SettableToPower<TDist>, CanGetAverageLog<TDist>,
								SettableToUniform, SettableTo<TDist>, SettableToRatio<TDist>, SettableToWeightedSum<TDist>
		{
#if DEBUG
			if (cases.Count != values.Count) throw new ArgumentException("cases.Count != values.Count");
#endif
			TDist uniform = (TDist)result.Clone();
			uniform.SetToUniform();
			return ExitAverageConditional2<TDist>(uniform, cases, values, result);
		}
		[Skip]
		public static TDist ExitAverageLogarithmInit<TDist>([IgnoreDependency] IList<TDist> values)
			where TDist : ICloneable
		{
			return (TDist)values[0].Clone();
		}
#else
		[Skip]
		public static DistributionArray<Bernoulli> CasesAverageLogarithm(DistributionArray<Bernoulli> result)
		{
			return result;
		}
		// result = prod_i values[i]^cases[i]  (messages out of a gate are blurred)
		public static T ExitAverageLogarithm<T>(DistributionArray<Bernoulli> cases, [SkipIfUniform] DistributionArray<T> values, T result)
			where T : Diffable, SettableTo<T>, ICloneable, SettableToUniform, SettableToProduct<T>,
			SettableToPower<T>, SettableToRatio<T>, SettableToWeightedSum<T>, LogInnerProductable<T>, CanGetAverageLog<T>
		{
			if (cases.Count != values.Count) throw new ArgumentException("cases.Count != values.Count");
			if (cases.Count == 0) throw new ArgumentException("cases.Count == 0");
			else {
				result.SetToPower(values[0], cases[0].LogOdds);
				if (cases.Count > 1) {
					// TODO: use pre-allocated buffer
					T power = (T)result.Clone();
					for (int i = 1; i < cases.Count; i++) {
						power.SetToPower(values[i], cases[i].LogOdds);
						result.SetToProduct(result, power);
					}
				}
			}
			return result;
		}
#endif
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.ExitTwo{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Gate), "ExitTwo<>")]
	[Quality(QualityBand.Experimental)]
	public static class GateExitTwoOp
	{
		/// <summary>
		/// EP message to 'values'
		/// </summary>
		/// <param name="exitTwo">Incoming message from 'exitTwo'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'values' as the random arguments are varied.
		/// The formula is <c>proj[p(values) sum_(exitTwo) p(exitTwo) factor(exitTwo,case0,case1,values)]/p(values)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="exitTwo"/> is not a proper distribution</exception>
		public static TList ValuesAverageConditional<T, TList>([SkipIfUniform] T exitTwo, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(exitTwo);
			return result;
		}

		// must takes values as input to distinguish from the other overload.
		/// <summary>
		/// EP message to 'case0'
		/// </summary>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'case0' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'case0' as the random arguments are varied.
		/// The formula is <c>proj[p(case0) sum_(values) p(values) factor(exitTwo,case0,case1,values)]/p(case0)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		[Skip]
		public static Bernoulli Case0AverageConditional<T>([SkipIfAllUniform] IList<T> values)
		{
			return Bernoulli.Uniform();
		}

		/// <summary>
		/// EP message to 'case1'
		/// </summary>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'case1' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'case1' as the random arguments are varied.
		/// The formula is <c>proj[p(case1) sum_(values) p(values) factor(exitTwo,case0,case1,values)]/p(case1)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		[Skip]
		public static Bernoulli Case1AverageConditional<T>([SkipIfAllUniform] IList<T> values)
		{
			return Bernoulli.Uniform();
		}
		/// <summary>
		/// EP message to 'exitTwo'
		/// </summary>
		/// <param name="exitTwo">Incoming message from 'exitTwo'.</param>
		/// <param name="case0">Incoming message from 'case0'.</param>
		/// <param name="case1">Incoming message from 'case1'.</param>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'exitTwo' as the random arguments are varied.
		/// The formula is <c>proj[p(exitTwo) sum_(case0,case1,values) p(case0,case1,values) factor(exitTwo,case0,case1,values)]/p(exitTwo)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		public static T ExitTwoAverageConditional<T>(T exitTwo, Bernoulli case0, Bernoulli case1, [SkipIfAllUniform] IList<T> values, T result)
			where T : SettableTo<T>, ICloneable, SettableToProduct<T>,
								SettableToRatio<T>, SettableToWeightedSum<T>
		{
			result.SetToProduct(exitTwo, values[0]);
			double scale = case0.LogOdds;
			double resultScale = scale;
			// TODO: use pre-allocated buffer
			T product = (T)exitTwo.Clone();
			product.SetToProduct(exitTwo, values[1]);
			scale = case1.LogOdds;
			double shift = Math.Max(resultScale, scale);
			// avoid (-Infinity) - (-Infinity)
			if (Double.IsNegativeInfinity(shift)) {
				throw new AllZeroException();
			} else {
				result.SetToSum(Math.Exp(resultScale - shift), result, Math.Exp(scale - shift), product);
				resultScale = MMath.LogSumExp(resultScale, scale);
			}
			result.SetToRatio(result, exitTwo);
			return result;
		}
		// Constant values array ////////////////////////
		/// <summary>
		/// EP message to 'case0'
		/// </summary>
		/// <param name="exitTwo">Incoming message from 'exitTwo'.</param>
		/// <param name="values">Incoming message from 'values'.</param>
		/// <returns>The outgoing EP message to the 'case0' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'case0' as the random arguments are varied.
		/// The formula is <c>proj[p(case0) sum_(exitTwo,values) p(exitTwo,values) factor(exitTwo,case0,case1,values)]/p(case0)</c>.
		/// </para></remarks>
		public static Bernoulli Case0AverageConditional<T, DomainType>(T exitTwo, IList<DomainType> values)
			where T : CanGetLogProb<DomainType>
		{
			return Bernoulli.FromLogOdds(exitTwo.GetLogProb(values[0]));
		}

		/// <summary>
		/// EP message to 'case1'
		/// </summary>
		/// <param name="exitTwo">Incoming message from 'exitTwo'.</param>
		/// <param name="values">Incoming message from 'values'.</param>
		/// <returns>The outgoing EP message to the 'case1' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'case1' as the random arguments are varied.
		/// The formula is <c>proj[p(case1) sum_(exitTwo,values) p(exitTwo,values) factor(exitTwo,case0,case1,values)]/p(case1)</c>.
		/// </para></remarks>
		public static Bernoulli Case1AverageConditional<T, DomainType>(T exitTwo, IList<DomainType> values)
			where T : CanGetLogProb<DomainType>
		{
			return Bernoulli.FromLogOdds(exitTwo.GetLogProb(values[1]));
		}

		/// <summary>
		/// EP message to 'exitTwo'
		/// </summary>
		/// <param name="exitTwo">Incoming message from 'exitTwo'.</param>
		/// <param name="case0">Incoming message from 'case0'.</param>
		/// <param name="case1">Incoming message from 'case1'.</param>
		/// <param name="values">Incoming message from 'values'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'exitTwo' as the random arguments are varied.
		/// The formula is <c>proj[p(exitTwo) sum_(case0,case1,values) p(case0,case1,values) factor(exitTwo,case0,case1,values)]/p(exitTwo)</c>.
		/// </para></remarks>
		public static T ExitTwoAverageConditional<T, DomainType>(T exitTwo, Bernoulli case0, Bernoulli case1, IList<DomainType> values, T result)
			where T : ICloneable, HasPoint<DomainType>, SettableToWeightedSum<T>
		{
			result.Point = values[0];
			double scale = case0.LogOdds;
			double resultScale = scale;
			// TODO: overload SetToSum to accept constants.
			T product = (T)exitTwo.Clone();
			product.Point = values[1];
			scale = case1.LogOdds;
			double shift = Math.Max(resultScale, scale);
			// avoid (-Infinity) - (-Infinity)
			if (Double.IsNegativeInfinity(shift)) {
				throw new AllZeroException();
				// do nothing
			} else {
				result.SetToSum(Math.Exp(resultScale - shift), result, Math.Exp(scale - shift), product);
				resultScale = MMath.LogSumExp(resultScale, scale);
			}
			return result;
		}

		//-- VMP ------------------------------------------------------------------------------------------------
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="exitTwo">Incoming message from 'exitTwo'.</param>
		/// <param name="case0">Incoming message from 'case0'.</param>
		/// <param name="case1">Incoming message from 'case1'.</param>
		/// <param name="values">Incoming message from 'values'.</param>
		/// <param name="to_exitTwo">Outgoing message to 'exitTwo'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor<T>(T exitTwo, Bernoulli case0, Bernoulli case1, IList<T> values, [Fresh] T to_exitTwo)
			where T : ICloneable, SettableToProduct<T>,
								SettableToPower<T>, CanGetAverageLog<T>,
								SettableToUniform, SettableTo<T>, SettableToRatio<T>, SettableToWeightedSum<T>
		{
			// cancel the evidence message from the child variable's child factors
			//T to_exit = (T)exitTwo.Clone();
			//to_exit = ExitTwoAverageLogarithm(case0, case1, values, to_exit);
			return -to_exitTwo.GetAverageLog(exitTwo);
		}
		// result is always exit (messages into a gate are unchanged)
		/// <summary>
		/// VMP message to 'values'
		/// </summary>
		/// <param name="exitTwo">Incoming message from 'exitTwo'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'values' with 'exitTwo' integrated out.
		/// The formula is <c>sum_exitTwo p(exitTwo) factor(exitTwo,case0,case1,values)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="exitTwo"/> is not a proper distribution</exception>
		public static TList ValuesAverageLogarithm<T, TList>([SkipIfUniform, Trigger] T exitTwo, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(exitTwo);
			return result;
		}

		/// <summary>
		/// VMP message to 'case0'
		/// </summary>
		/// <param name="exitTwo">Incoming message from 'exitTwo'.</param>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'case0' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'case0'.
		/// Because the factor is deterministic, 'exitTwo' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(values) p(values) log(sum_exitTwo p(exitTwo) factor(exitTwo,case0,case1,values)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		public static Bernoulli Case0AverageLogarithm<T>(T exitTwo, [SkipIfAllUniform, Proper] IList<T> values)
			where T : CanGetAverageLog<T>
		{
			return Bernoulli.FromLogOdds(values[0].GetAverageLog(exitTwo));
		}

		/// <summary>
		/// VMP message to 'case1'
		/// </summary>
		/// <param name="exitTwo">Incoming message from 'exitTwo'.</param>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'case1' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'case1'.
		/// Because the factor is deterministic, 'exitTwo' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(values) p(values) log(sum_exitTwo p(exitTwo) factor(exitTwo,case0,case1,values)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		public static Bernoulli Case1AverageLogarithm<T>(T exitTwo, [SkipIfAllUniform, Proper] IList<T> values)
			where T : CanGetAverageLog<T>
		{
			return Bernoulli.FromLogOdds(values[1].GetAverageLog(exitTwo));
		}


		// result = prod_i values[i]^cases[i]  (messages out of a gate are blurred)
		/// <summary>
		/// VMP message to 'exitTwo'
		/// </summary>
		/// <param name="case0">Incoming message from 'case0'.</param>
		/// <param name="case1">Incoming message from 'case1'.</param>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'exitTwo' as the random arguments are varied.
		/// The formula is <c>proj[sum_(case0,case1,values) p(case0,case1,values) factor(exitTwo,case0,case1,values)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		public static T ExitTwoAverageLogarithm<T>(Bernoulli case0, Bernoulli case1, [SkipIfAllUniform, Proper] IList<T> values, T result)
			where T : ICloneable, SettableToProduct<T>,
								SettableToPower<T>, CanGetAverageLog<T>,
								SettableToUniform, SettableTo<T>, SettableToRatio<T>, SettableToWeightedSum<T>
		{
			T uniform = (T)result.Clone();
			uniform.SetToUniform();
			return ExitTwoAverageConditional<T>(uniform, case0, case1, values, result);
		}
	}


	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.ExitRandom{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Gate), "ExitRandom<>")]
	[Quality(QualityBand.Preview)]
	public static class GateExitRandomOp
	{
		[Skip]
		public static double LogAverageFactor()		{			return 0.0;		}

		/// <summary>
		/// Gibbs message to 'values'.
		/// </summary>
		/// <param name="exit">Incoming point message from 'Exit'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		public static T[] ValuesAverageConditional<T>([IsReturnedInEveryElement] T exit, T[] result)
		{
			for (int i=0; i < result.Length; i++)
				result[i] = exit;
			return result;
		}

		// result is always exit (messages into a gate are unchanged)
		/// <summary>
		/// Gibbs message to 'values'.
		/// </summary>
		/// <param name="exit">Incoming message from 'Exit'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'values'.
		/// The formula is <c>int log(f(values,x)) q(x) dx</c> where <c>x = (Exit,cases)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="exit"/> is not a proper distribution</exception>
		public static TList ValuesAverageConditional<T, TList>([IsReturnedInEveryElement] T exit, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			return ValuesAverageLogarithm(exit, result);
		}
		/// <summary>
		/// Gibbs message to 'cases'.
		/// </summary>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		[Skip]
		public static BernoulliList CasesAverageConditional<BernoulliList>(BernoulliList result)
			where BernoulliList : SettableToUniform
		{
			return CasesAverageLogarithm<BernoulliList>(result);
		}

		/// <summary>
		/// Gibbs message to 'Exit'
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="cases"></param>
		/// <param name="values"></param>
		/// <returns></returns>
		public static T ExitAverageConditional<T>(bool[] cases, [SkipIfAllUniform] IList<T> values)
		{
			for (int i = 0; i < cases.Length; i++) {
				if (cases[i]) return values[i];
			}
			// cases was entirely false
			throw new ArgumentException("cases is all false");
		}

		///// <summary>
		///// Gibbs message to 'Exit'
		///// </summary>
		///// <typeparam name="T"></typeparam>
		///// <typeparam name="TDomain"></typeparam>
		///// <param name="cases"></param>
		///// <param name="values"></param>
		///// <param name="result"></param>
		///// <returns></returns>
		//public static T ExitAverageConditional<T, TDomain>(IList<bool> cases, IList<TDomain> values, T result)
		//    where T : IDistribution<TDomain>
		//{
		//    for (int i = 0; i < cases.Count; i++)
		//    {
		//        if (cases[i]) return Distribution.SetPoint<T,TDomain>(result, values[i]);
		//    }
		//    // cases was entirely false
		//    throw new ArgumentException("cases is all false");
		//}

		//-- VMP ------------------------------------------------------------------------------------------------
		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (Exit,cases,values)</c>.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor()
		{
			return 0.0;
		}
		// result is always exit (messages into a gate are unchanged)
		/// <summary>
		/// VMP message to 'values'.
		/// </summary>
		/// <param name="exit">Incoming message from 'Exit'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'values'.
		/// The formula is <c>int log(f(values,x)) q(x) dx</c> where <c>x = (Exit,cases)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="exit"/> is not a proper distribution</exception>
		public static TList ValuesAverageLogarithm<T, TList>([IsReturnedInEveryElement] T exit, TList result)
			where TList : CanSetAllElementsTo<T>
		{
			result.SetAllElementsTo(exit);
			return result;
		}
		/// <summary>
		/// VMP message to 'cases'.
		/// </summary>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'cases'.
		/// The formula is <c>int log(f(cases,x)) q(x) dx</c> where <c>x = (Exit,values)</c>.
		/// </para></remarks>
		[Skip]
		public static BernoulliList CasesAverageLogarithm<BernoulliList>(BernoulliList result)
			where BernoulliList : SettableToUniform
		{
			result.SetToUniform();
			return result;
		}

		// result = prod_i values[i]^cases[i]  (messages out of a gate are blurred)
		/// <summary>
		/// VMP message to 'Exit'.
		/// </summary>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="values">Incoming message from 'values'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'Exit'.
		/// The formula is <c>int log(f(Exit,x)) q(x) dx</c> where <c>x = (cases,values)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="values"/> is not a proper distribution</exception>
		public static T ExitAverageLogarithm<T>(IList<Bernoulli> cases, [SkipIfAllUniform, Proper] IList<T> values, T result)
			where T : ICloneable, SettableToProduct<T>,
			SettableToPower<T>, CanGetAverageLog<T>,
			SettableToUniform, SettableTo<T>, SettableToRatio<T>, SettableToWeightedSum<T>
		{
			if (cases.Count != values.Count) throw new ArgumentException("cases.Count != values.Count");
			if (cases.Count == 0) throw new ArgumentException("cases.Count == 0");
			else {
				result.SetToPower(values[0], Math.Exp(cases[0].LogOdds));
				if (cases.Count > 1) {
					// TODO: use pre-allocated buffer
					T power = (T)result.Clone();
					for (int i = 1; i < cases.Count; i++) {
						power.SetToPower(values[i], Math.Exp(cases[i].LogOdds));
						result.SetToProduct(result, power);
					}
				}
			}
			return result;
		}
		[Skip]
		public static T ExitAverageLogarithmInit<T>([IgnoreDependency] IList<T> values)
			where T : ICloneable
		{
			return (T)values[0].Clone();
		}
	}
}
