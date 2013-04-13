// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Factors for handling gates.
	/// </summary>
	[Hidden]
	public static class Gate
	{
		/// <summary>
		/// Enter factor
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="selector"></param>
		/// <param name="value"></param>
		/// <returns></returns>
		public static T[] Enter<T>(int selector, [IsReturnedInEveryElement] T value)
		{
			throw new NotImplementedException();
			T[] result = new T[2];
			for (int i = 0; i < 2; i++) {
				result[i] = value;
			}
			return result;
		}
		public static T[] Enter<T>(bool selector, [IsReturnedInEveryElement] T value)
		{
			T[] result = new T[2];
			for (int i = 0; i < 2; i++) {
				result[i] = value;
			}
			return result;
		}

		/// <summary>
		/// Enter partial factor
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="selector"></param>
		/// <param name="value"></param>
		/// <param name="indices"></param>
		/// <returns></returns>
		public static T[] EnterPartial<T>(int selector, [IsReturnedInEveryElement] T value, int[] indices)
		{
			T[] result = new T[indices.Length];
			for (int i = 0; i < indices.Length; i++) {
				result[i] = value;
			}
			return result;
		}
		public static T[] EnterPartial<T>(bool selector, [IsReturnedInEveryElement] T value, int[] indices)
		{
			T[] result = new T[indices.Length];
			for (int i = 0; i < indices.Length; i++) {
				result[i] = value;
			}
			return result;
		}

		/// <summary>
		/// Enter partial factor with two cases
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="case0"></param>
		/// <param name="case1"></param>
		/// <param name="value"></param>
		/// <param name="indices"></param>
		/// <returns></returns>
		public static T[] EnterPartialTwo<T>(bool case0, bool case1, [IsReturnedInEveryElement] T value, int[] indices)
		{
			T[] result = new T[indices.Length];
			for (int i = 0; i < indices.Length; i++)
			{
				result[i] = value;
			}
			return result;
		}

		/// <summary>
		/// Enter one factor
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="selector"></param>
		/// <param name="value"></param>
		/// <param name="index"
		/// <returns></returns>
		public static T EnterOne<T>(int selector, [IsReturned] T value, int index)
		{
			return value;
		}

		/// <summary>
		/// Exit factor
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="cases"></param>
		/// <param name="values"></param>
		/// <returns></returns>
		public static T Exit<T>(bool[] cases, T[] values)
		{
			for (int i = 0; i < cases.Length; i++)
				if (cases[i])
					return values[i];

			throw new ApplicationException("Exit factor: no case is true");
		}

		/// <summary>
		/// Exit factor with two cases
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="case0"></param>
		/// <param name="case1"></param>
		/// <param name="values"></param>
		/// <returns></returns>
		public static T ExitTwo<T>(bool case0, bool case1, T[] values)
		{
			if (case0)
				return values[0];
			else if (case1)
				return values[1];

			throw new ApplicationException("ExitTwo factor: neither case is true");
		}

		/// <summary>
		/// Exit random factor
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="cases"></param>
		/// <param name="values"></param>
		/// <returns></returns>
		[Stochastic]
		[ParameterNames("Exit", "cases", "values")]
		public static T ExitRandom<T>(bool[] cases, T[] values)
		{
			return Exit(cases, values);
		}

		/// <summary>
		/// Exit random factor
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="case0"></param>
		/// <param name="case1"></param>
		/// <param name="values"></param>
		/// <returns></returns>
		[Stochastic]
		[ParameterNames("Exit", "cases", "values")]
		public static T ExitRandomTwo<T>(bool case0, bool case1, T[] values)
		{
			if (case0)
				return values[0];
			else if (case1)
				return values[1];

			throw new ApplicationException("ExitTwo factor: neither case is true");
		}

#if true
		/// <summary>
		/// Exiting variable
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="Def"></param>
		/// <param name="Marginal"></param>
		/// <returns></returns>
		[ParameterNames("Use", "Def", "Marginal")]
		public static T ExitingVariable<T>(T Def, out T Marginal)
		{
			throw new InvalidOperationException("Should never be called with deterministic arguments");
		}
		/// <summary>
		/// Replicate an exiting variable
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="Def"></param>
		/// <returns></returns>
		[ParameterNames("Uses", "Def", "count")]
		public static T[] ReplicateExiting<T>(T Def, int count)
		{
			throw new InvalidOperationException("Should never be called with deterministic arguments");
		}
#else
		/// <summary>
		/// Exiting variable factor
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="Def"></param>
		/// <param name="Marginal"></param>
		/// <returns></returns>
		[ParameterNames("Uses", "Def", "Marginal")]
		public static T[] ExitingVariable<T>(T Def, T Marginal)
		{
			throw new InvalidOperationException("Should never be called with deterministic arguments");
		}
#endif

		/// <summary>
		/// Boolean cases factor
		/// </summary>
		/// <param name="b"></param>
		/// <returns></returns>
		public static bool[] Cases(bool b)
		{
			bool[] result = new bool[2];
			result[0] = b;
			result[1] = !b;
			return result;
		}

		/// <summary>
		/// Boolean cases factor expanded into elements
		/// </summary>
		/// <param name="b"></param>
		/// <param name="case0">case 0 (true)</param>
		/// <param name="case1">case 1 (false)</param>
		/// <returns></returns>
		public static void CasesBool(bool b, out bool case0, out bool case1)
		{
			case0 = b;
			case1 = !b;
		}

		// TODO: fix bug which prevents this being called 'Cases'
		/// <summary>
		/// Integer cases factor
		/// </summary>
		/// <param name="i">index</param>
		/// <param name="count">number of cases</param>
		/// <returns></returns>
		public static bool[] CasesInt(int i, int count)
		{
			bool[] result = new bool[count];
			for (int j = 0; j < count; j++)
				result[j] = false;
			result[i] = true;
			return result;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.Cases"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Gate), "Cases", typeof(bool))]
	[Quality(QualityBand.Mature)]
	public static class CasesOp
	{
		// result.LogOdds = [log p(b=true), log p(b=false)]
		/// <summary>
		/// EP message to 'cases'.
		/// </summary>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'cases'.
		/// The formula is <c>int f(cases,x) q(x) dx</c> where <c>x = (b)</c>.
		/// </para></remarks>
		public static BernoulliList CasesAverageConditional<BernoulliList>(Bernoulli b, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
		{
			if (result.Count != 2) throw new ArgumentException("result.Count != 2");
			result[0] = Bernoulli.FromLogOdds(b.GetLogProbTrue());
			result[1] = Bernoulli.FromLogOdds(b.GetLogProbFalse());
			return result;
		}
		[Skip]
		public static DistributionStructArray<Bernoulli, bool> CasesAverageConditionalInit()
		{
			return new DistributionStructArray<Bernoulli, bool>(2);
		}

		// result = p(b=true) / (p(b=true) + p(b=false))
		//        = 1 / (1 + p(b=false)/p(b=true))
		//        = 1 / (1 + exp(-(log p(b=true) - log p(b=false)))
		// where cases[0].LogOdds = log p(b=true)
		//       cases[1].LogOdds = log p(b=false)
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="cases">Incoming message from 'cases'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (cases)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="cases"/> is not a proper distribution</exception>
		public static Bernoulli BAverageConditional([SkipIfAnyUniform] IList<Bernoulli> cases)
		{
			// avoid (-Infinity) - (-Infinity)
			if (cases[0].LogOdds == cases[1].LogOdds)
			{
				if (Double.IsNegativeInfinity(cases[0].LogOdds) && Double.IsNegativeInfinity(cases[1].LogOdds)) throw new AllZeroException();
				return new Bernoulli();
			}
			else
			{
				return Bernoulli.FromLogOdds(cases[0].LogOdds - cases[1].LogOdds);
			}
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="cases">Incoming message from 'cases'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		public static double LogEvidenceRatio(IList<Bernoulli> cases, Bernoulli b)
		{
			// result = log (p(data|b=true) p(b=true) + p(data|b=false) p(b=false))
			//          log (p(data|b=true) p(b=true) + p(data|b=false) (1-p(b=true))
			//          log ((p(data|b=true) - p(data|b=false)) p(b=true) + p(data|b=false))
			//          log ((p(data|b=true)/p(data|b=false) - 1) p(b=true) + 1) + log p(data|b=false)
			// where cases[0].LogOdds = log p(data|b=true)
			//       cases[1].LogOdds = log p(data|b=false)
			if (b.IsPointMass) return b.Point ? cases[0].LogOdds : cases[1].LogOdds;
			//else return MMath.LogSumExp(cases[0].LogOdds + b.GetLogProbTrue(), cases[1].LogOdds + b.GetLogProbFalse());
			else {
				// the common case is when cases[0].LogOdds == cases[1].LogOdds.  we must not introduce rounding error in that case.
				if (cases[0].LogOdds >= cases[1].LogOdds) {
					if (Double.IsNegativeInfinity(cases[1].LogOdds))
						return cases[0].LogOdds + b.GetLogProbTrue();
					else
						return cases[1].LogOdds + MMath.Log1Plus(b.GetProbTrue() * MMath.ExpMinus1(cases[0].LogOdds - cases[1].LogOdds));
				} else {
					if (Double.IsNegativeInfinity(cases[0].LogOdds))
						return cases[1].LogOdds + b.GetLogProbFalse();
					else
						return cases[0].LogOdds + MMath.Log1Plus(b.GetProbFalse() * MMath.ExpMinus1(cases[1].LogOdds - cases[0].LogOdds));
				}
			}
		}

		// Gibbs message
		[Skip]
		public static double LogEvidenceRatio(IList<Bernoulli> cases, bool b)
		{
			return 0.0;
			//return b ? cases[0].LogOdds : cases[1].LogOdds;
		}

		//-- VMP --------------------------------------------------------------------------------------------

		/// <summary>
		/// VMP message to 'cases'.
		/// </summary>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'cases'.
		/// The formula is <c>int log(f(cases,x)) q(x) dx</c> where <c>x = (b)</c>.
		/// </para></remarks>
		public static BernoulliList CasesAverageLogarithm<BernoulliList>(Bernoulli b, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
		{
			return CasesAverageConditional(b, result);
		}
		[Skip]
		public static DistributionStructArray<Bernoulli, bool> CasesAverageLogarithmInit()
		{
			return new DistributionStructArray<Bernoulli, bool>(2);
		}

		[Skip]
		public static DistributionType CasesDeriv<DistributionType>(DistributionType result)
			where DistributionType : SettableToUniform
		{
			result.SetToUniform();
			return result;
		}

		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="cases">Incoming message from 'cases'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (cases)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="cases"/> is not a proper distribution</exception>
		// TM: SkipIfAny (rather than SkipIfAll) is important for getting good schedules
		public static Bernoulli BAverageLogarithm([SkipIfAnyUniform] IList<Bernoulli> cases)
		{
			return BAverageConditional(cases);
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="cases">Incoming message from 'cases'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (cases,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="cases"/> is not a proper distribution</exception>
		public static double AverageLogFactor([SkipIfAnyUniform] IList<Bernoulli> cases, Bernoulli b)
		{
			double probTrue = b.GetProbTrue();
			return probTrue * cases[0].LogOdds + (1 - probTrue) * cases[1].LogOdds;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.CasesBool"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Gate), "CasesBool", typeof(bool), typeof(bool), typeof(bool))]
	[Quality(QualityBand.Experimental)]
	public static class CasesBoolOp
	{
		/// <summary>
		/// EP message to 'case0'
		/// </summary>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'case0' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'case0' as the random arguments are varied.
		/// The formula is <c>proj[p(case0) sum_(b) p(b) factor(b,case0,case1)]/p(case0)</c>.
		/// </para></remarks>
		public static Bernoulli Case0AverageConditional(Bernoulli b)
		{
			return Bernoulli.FromLogOdds(b.GetLogProbTrue());
		}

		/// <summary>
		/// EP message to 'case1'
		/// </summary>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'case1' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'case1' as the random arguments are varied.
		/// The formula is <c>proj[p(case1) sum_(b) p(b) factor(b,case0,case1)]/p(case1)</c>.
		/// </para></remarks>
		public static Bernoulli Case1AverageConditional(Bernoulli b)
		{
			return Bernoulli.FromLogOdds(b.GetLogProbFalse());
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="case0">Incoming message from 'case0'.</param>
		/// <param name="case1">Incoming message from 'case1'.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(case0,case1) p(case0,case1) factor(b,case0,case1)]/p(b)</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		public static Bernoulli BAverageConditional(Bernoulli case0, Bernoulli case1)
		{
			// result = p(b=true) / (p(b=true) + p(b=false))
			//        = 1 / (1 + p(b=false)/p(b=true))
			//        = 1 / (1 + exp(-(log p(b=true) - log p(b=false)))
			// where cases[0].LogOdds = log p(b=true)
			//       cases[1].LogOdds = log p(b=false)
			// avoid (-Infinity) - (-Infinity)
			if (case0.LogOdds == case1.LogOdds)
			{
				if (Double.IsNegativeInfinity(case0.LogOdds) && Double.IsNegativeInfinity(case1.LogOdds)) throw new AllZeroException();
				return new Bernoulli();
			}
			else
			{
				return Bernoulli.FromLogOdds(case0.LogOdds - case1.LogOdds);
			}
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="case0">Incoming message from 'case0'.</param>
		/// <param name="case1">Incoming message from 'case1'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(case0,case1,b) p(case0,case1,b) factor(b,case0,case1))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Bernoulli case0, Bernoulli case1, Bernoulli b)
		{
			// result = log (p(data|b=true) p(b=true) + p(data|b=false) p(b=false))
			// where cases[0].LogOdds = log p(data|b=true)
			//       cases[1].LogOdds = log p(data|b=false)
			if (b.IsPointMass) return b.Point ? case0.LogOdds : case1.LogOdds;
			else return MMath.LogSumExp(case0.LogOdds + b.GetLogProbTrue(), case1.LogOdds + b.GetLogProbFalse());
		}

		//-- VMP --------------------------------------------------------------------------------------------

		/// <summary>
		/// VMP message to 'case0'
		/// </summary>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'case0' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'case0'.
		/// The formula is <c>exp(sum_(b) p(b) log(factor(b,case0,case1)))</c>.
		/// </para></remarks>
		public static Bernoulli Case0AverageLogarithm(Bernoulli b)
		{
			return Case0AverageConditional(b);
		}

		/// <summary>
		/// VMP message to 'case1'
		/// </summary>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'case1' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'case1'.
		/// The formula is <c>exp(sum_(b) p(b) log(factor(b,case0,case1)))</c>.
		/// </para></remarks>
		public static Bernoulli Case1AverageLogarithm(Bernoulli b)
		{
			return Case1AverageConditional(b);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="case0">Incoming message from 'case0'.</param>
		/// <param name="case1">Incoming message from 'case1'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// The formula is <c>exp(sum_(case0,case1) p(case0,case1) log(factor(b,case0,case1)))</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		public static Bernoulli BAverageLogarithm(Bernoulli case0, Bernoulli case1)
		{
			return BAverageConditional(case0, case1);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="case0">Incoming message from 'case0'.</param>
		/// <param name="case1">Incoming message from 'case1'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[SkipIfAllUniform("case0", "case1")]
		public static double AverageLogFactor(Bernoulli case0, Bernoulli case1, Bernoulli b)
		{
			double probTrue = b.GetProbTrue();
			return probTrue * case0.LogOdds + (1 - probTrue) * case1.LogOdds;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Gate.CasesInt"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(new string[] { "Cases", "i", "count" }, typeof(Gate), "CasesInt", typeof(int), typeof(int))]
	[Quality(QualityBand.Mature)]
	public static class IntCasesOp
	{
		/// <summary>
		/// EP message to 'casesInt'.
		/// </summary>
		/// <param name="i">Incoming message from 'i'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'casesInt'.
		/// The formula is <c>int f(casesInt,x) q(x) dx</c> where <c>x = (i)</c>.
		/// </para></remarks>
		public static BernoulliList CasesAverageConditional<BernoulliList>(Discrete i, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
		{
			if (result.Count != i.Dimension) throw new ArgumentException("result.Count ("+result.Count+") != i.Dimension ("+i.Dimension+")");
			for (int j = 0; j < result.Count; j++)
			{
				result[j] = Bernoulli.FromLogOdds(i.GetLogProb(j));
			}
			return result;
		}
		[Skip]
		public static DistributionStructArray<Bernoulli, bool> CasesAverageConditionalInit([IgnoreDependency] Discrete i)
		{
			return new DistributionStructArray<Bernoulli, bool>(i.Dimension);
		}

		/// <summary>
		/// EP message to 'i'.
		/// </summary>
		/// <param name="cases">Incoming message from 'casesInt'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'i'.
		/// The formula is <c>int f(i,x) q(x) dx</c> where <c>x = (casesInt)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="cases"/> is not a proper distribution</exception>
		public static Discrete IAverageConditional([SkipIfUniform] IList<Bernoulli> cases, Discrete result)
		{
			if (cases.Count != result.Dimension) throw new ArgumentException("cases.Count ("+cases.Count+") != result.Dimension ("+result.Dimension+")");
			Vector probs = result.GetWorkspace();
			double max = cases[0].LogOdds;
			for (int j = 1; j < cases.Count; j++)
			{
				if (cases[j].LogOdds > max) max = cases[j].LogOdds;
			}
			// avoid (-Infinity) - (-Infinity)
			if (Double.IsNegativeInfinity(max)) throw new AllZeroException();
            // TODO: how to avoid having this specific line
            if (probs.Sparsity.IsApproximate) probs.SetAllElementsTo(0);
			for (int j = 0; j < result.Dimension; j++)
			{
				probs[j] = Math.Exp(cases[j].LogOdds - max);
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="cases">Incoming message from 'casesInt'.</param>
		public static double LogEvidenceRatio([SkipIfUniform] IList<Bernoulli> cases, Discrete i)
		{
			if (i.IsPointMass) return cases[i.Point].LogOdds;
			else
			{
				double[] logOdds = new double[cases.Count];
				for (int j = 0; j < cases.Count; j++)
				{
					logOdds[j] = cases[j].LogOdds + i.GetLogProb(j);
				}
				return MMath.LogSumExp(logOdds);
			}
		}
		/// <summary>
		/// Evidence message for Gibbs.
		/// </summary>
		/// <param name="cases">Incoming message from 'casesInt'.</param>
		public static double LogEvidenceRatio(IList<Bernoulli> cases, int i)
		{
			return cases[i].LogOdds;
		}

		//-- VMP --------------------------------------------------------------------------------------------
		/// <summary>
		/// VMP message to 'casesInt'.
		/// </summary>
		/// <param name="i">Incoming message from 'i'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'casesInt'.
		/// The formula is <c>int log(f(casesInt,x)) q(x) dx</c> where <c>x = (i)</c>.
		/// </para></remarks>
		public static BernoulliList CasesAverageLogarithm<BernoulliList>(Discrete i, BernoulliList result)
			where BernoulliList : IList<Bernoulli>
		{
			return CasesAverageConditional(i, result);
		}
		[Skip]
		public static DistributionStructArray<Bernoulli, bool> CasesAverageLogarithmInit([IgnoreDependency] Discrete i)
		{
			return new DistributionStructArray<Bernoulli, bool>(i.Dimension);
		}
		/// <summary>
		/// VMP message to 'i'.
		/// </summary>
		/// <param name="cases">Incoming message from 'casesInt'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'i'.
		/// The formula is <c>int log(f(i,x)) q(x) dx</c> where <c>x = (casesInt)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="cases"/> is not a proper distribution</exception>
		public static Discrete IAverageLogarithm([SkipIfUniform] IList<Bernoulli> cases, Discrete result)
		{
			return IAverageConditional(cases, result);
		}
		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="cases">Incoming message from 'casesInt'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="i">Incoming message from 'i'.</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (casesInt,i)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="cases"/> is not a proper distribution</exception>
		public static double AverageLogFactor([SkipIfAllUniform] IList<Bernoulli> cases, Discrete i)
		{
			double sum = 0.0;
			for (int j = 0; j < cases.Count; j++)
			{
				sum += i[j] * cases[j].LogOdds;
			}
			return sum;
		}
	}
}
