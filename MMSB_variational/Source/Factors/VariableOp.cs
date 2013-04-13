// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing EP messages for <see cref="Factor.Variable&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Variable<>")]
	[FactorMethod(typeof(Factor), "VariableInit<>")]
	[Quality(QualityBand.Preview)]
	public static class VariableOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(use,def,marginal))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio()
		{
			return 0.0;
		}

		/// <summary>
		/// EP message to 'marginal'
		/// </summary>
		/// <param name="use">Incoming message from 'use'.</param>
		/// <param name="def">Incoming message from 'def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'marginal' as the random arguments are varied.
		/// The formula is <c>proj[p(marginal) sum_(use,def) p(use,def) factor(use,def,marginal)]/p(marginal)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="def"/> is not a proper distribution</exception>
		[SkipIfAllUniform]
		[MultiplyAll]
		public static T MarginalAverageConditional<T>(T use, T def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			result.SetToProduct(def, use);
			return result;
		}
		[Skip]
		public static T MarginalAverageConditionalInit<T>([IgnoreDependency] T def)
			where T: ICloneable
		{
			return (T)def.Clone();
		}

		/// <summary>
		/// EP message to 'use'
		/// </summary>
		/// <param name="Def">Incoming message from 'def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'use' as the random arguments are varied.
		/// The formula is <c>proj[p(use) sum_(def) p(def) factor(use,def,marginal)]/p(use)</c>.
		/// </para></remarks>
		public static T UseAverageConditional<T>([IsReturned] T Def)
		{
			return Def;
		}

		/// <summary>
		/// EP message to 'def'
		/// </summary>
		/// <param name="use">Incoming message from 'use'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'def' as the random arguments are varied.
		/// The formula is <c>proj[p(def) sum_(use) p(use) factor(use,def,marginal)]/p(def)</c>.
		/// </para></remarks>
		public static T DefAverageConditional<T>([IsReturned] T use)
		{
			return use;
		}
	}
	[FactorMethod(typeof(Factor), "VariableInit<>", Default=true)]
	[Quality(QualityBand.Preview)]
	public static class VariableInitOp
	{
	}

	/// <summary>
	/// Provides outgoing Gibbs messages for <see cref="Factor.Variable&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "VariableGibbs<>")]
	[Quality(QualityBand.Preview)]
	public static class VariableGibbsOp
	{
		/// <summary>
		/// Gibbs evidence
		/// </summary>
		/// <returns></returns>
		public static double GibbsEvidence<TDist,T>(TDist Use, TDist Def, GibbsMarginal<TDist, T> marginal)
			where TDist : IDistribution<T>, Sampleable<T>, CanGetLogAverageOf<TDist>
		{
			return Def.GetLogAverageOf(Use) -Def.GetLogProb(marginal.LastSample) -Use.GetLogProb(marginal.LastSample);
		}

		/// <summary>
		/// Gibbs message to 'Marginal'.
		/// </summary>
		/// <param name="Use">Incoming message from 'use'.</param>
		/// <param name="Def">Incoming message from 'def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the product of 'Def' and 'Uses' messages.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		[Stochastic]
		public static GibbsMarginal<TDist,T> MarginalGibbs<TDist, T>(
			TDist Use,
			[SkipIfUniform] TDist Def,
			GibbsMarginal<TDist, T> to_marginal) // must not be called 'result', because its value is used
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableTo<TDist>, Sampleable<T>
		{
			GibbsMarginal<TDist,T> result = to_marginal;
			TDist marginal = result.LastConditional;
			marginal.SetToProduct(Def, Use);
			result.LastConditional = marginal;
			// Allow a sample to be drawn from the last conditional, and add it to the sample
			// list and conditional list
			result.PostUpdate();
			return result;
		}
		[Stochastic]
		public static GibbsMarginal<TDist, T> MarginalGibbs<TDist, T>(
			T Use,
			[SkipIfUniform] TDist Def,
			GibbsMarginal<TDist, T> to_marginal) // must not be called 'result', because its value is used
			where TDist : IDistribution<T>, Sampleable<T>
		{
			GibbsMarginal<TDist,T> result = to_marginal;
			TDist marginal = result.LastConditional;
			marginal.Point = Use;
			result.LastConditional = marginal;
			// Allow a sample to be drawn from the last conditional, and add it to the sample
			// list and conditional list
			result.PostUpdate();
			return result;
		}

		/// <summary>
		/// Gibbs sample message to 'Use'
		/// </summary>
		/// <param name="marginal">Incoming message from 'marginal'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the current Gibbs sample.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="marginal"/> is not a proper distribution</exception>
		public static T UseGibbs<TDist, T>([SkipIfUniform] GibbsMarginal<TDist, T> marginal, T result)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			return marginal.LastSample;
		}

		/// <summary>
		/// Gibbs distribution message to 'Def'
		/// </summary>
		/// <param name="Def">Incoming message from 'def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the product of the 'Def' message with all 'Uses' messages except the current
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		public static TDist UseGibbs<TDist, T>([IsReturned] TDist Def, TDist result)
			where TDist : SettableTo<TDist>
		{
			result.SetTo(Def);
			return result;
		}

		/// <summary>
		/// Gibbs sample message to 'Def'
		/// </summary>
		/// <param name="marginal">Incoming message from 'marginal'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the current Gibbs sample.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="marginal"/> is not a proper distribution</exception>
		public static T DefGibbs<TDist, T>([SkipIfUniform] GibbsMarginal<TDist, T> marginal, T result)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			return marginal.LastSample;
		}
		public static TDist DefGibbs<TDist, T>([IsReturned] TDist Use, TDist result)
			where TDist : IDistribution<T>, SettableTo<TDist>
		{
			result.SetTo(Use);
			return result;
		}
	}

	/// <summary>
	/// Provides outgoing max product messages for <see cref="Factor.Variable&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "VariableMax<>")]
	[Quality(QualityBand.Preview)]
	public static class VariableMaxOp
	{
		public static T UseMaxConditional<T>(T Def, T result)
			where T : SettableTo<T>
		{
			result.SetTo(Def);
			if (result is UnnormalizedDiscrete) ((UnnormalizedDiscrete)(object)result).SetMaxToZero();
			return result;
		}
		[Skip]
		public static T UseMaxConditionalInit<T>([IgnoreDependency] T Def)
			where T : ICloneable
		{
			return (T)Def.Clone();
		}
		public static T DefMaxConditional<T>([IsReturned] T Use)
		{
			return Use;
		}
		[MultiplyAll]
		public static T MarginalMaxConditional<T>(T Use, [SkipIfUniform] T Def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			T res = VariableOp.MarginalAverageConditional<T>(Use, Def, result);
			if (res is UnnormalizedDiscrete) ((UnnormalizedDiscrete)(object)res).SetMaxToZero();
			return res;
		}
		[Skip]
		public static T MarginalMaxConditionalInit<T>([IgnoreDependency] T def)
			where T : ICloneable
		{
			return (T)def.Clone();
		}
	}

#if true
	/// <summary>
	/// Provides outgoing VMP messages for <see cref="Factor.Variable&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Variable<>")]
	[FactorMethod(typeof(Factor), "VariableInit<>")]
	[Quality(QualityBand.Preview)]
	public static class VariableVmpOp
	{
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="marginal">Outgoing message to 'marginal'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Uses,Def,Marginal))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor<T>([SkipIfUniform] T to_marginal /*, [IgnoreDependency] T def*/)
			where T : CanGetAverageLog<T>
		{
			return -to_marginal.GetAverageLog(to_marginal);
		}

		[SkipIfAllUniform]
		[MultiplyAll]
		public static T MarginalAverageLogarithm<T>(T use, [Proper] T def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			result.SetToProduct(def, use);
			return result;
		}
		[Skip]
		public static T MarginalAverageLogarithmInit<T>([IgnoreDependency] T def)
				where T : ICloneable
		{
			return (T)def.Clone();
		}

		/// <summary>
		/// VMP message to 'use'
		/// </summary>
		/// <param name="marginal">Current 'marginal'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'use' conditioned on the given values.
		/// </para></remarks>
		public static T UseAverageLogarithm<T>([IsReturned] T to_marginal, T result)
			where T : SettableTo<T>
		{
			result.SetTo(to_marginal);
			return result;
		}

		[Skip]
		public static T UseAverageLogarithmInit<T>([IgnoreDependency] T def)
			where T: ICloneable
		{
			return (T)def.Clone();
		}

		/// <summary>
		/// VMP message to 'def'
		/// </summary>
		/// <param name="marginal">Current 'marginal'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'def' conditioned on the given values.
		/// </para></remarks>
		public static T DefAverageLogarithm<T>([IsReturned] T to_marginal, T result)
			where T : SettableTo<T>
		{
			result.SetTo(to_marginal);
			return result;
		}
	}
#else
	/// <summary>
	/// Provides outgoing VMP messages for <see cref="Factor.Variable&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Variable<>")]
	[FactorMethod(typeof(Factor), "VariableInit<>", Default=false)]
	[Buffers("marginalB")]
	[Quality(QualityBand.Preview)]
	public static class VariableVmpBufferOp
	{
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="marginal">Outgoing message to 'marginal'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Uses,Def,Marginal))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor<T>([Fresh, SkipIfUniform] T marginalB, [IgnoreDependency] T def)
			where T : CanGetAverageLog<T>
		{
			return -marginalB.GetAverageLog(marginalB);
		}

		[SkipIfAllUniform]
		[MultiplyAll]
		public static T MarginalB<T>(T use, T def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			result.SetToProduct(def, use);
			return result;
		}
		//[Skip]
		//public static T MarginalBInit<T>([IgnoreDependency] T def)
		//  where T : ICloneable
		//{
		//  return (T)def.Clone();
		//}

		public static T MarginalAverageLogarithm<T>([IsReturned] T marginalB, T result)
			where T : SettableTo<T>
		{
			result.SetTo(marginalB);
			return result;
		}

		/// <summary>
		/// VMP message to 'use'
		/// </summary>
		/// <param name="marginal">Current 'marginal'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'use' conditioned on the given values.
		/// </para></remarks>
		public static T UseAverageLogarithm<T>([IsReturned] T marginalB, T result)
			where T : SettableTo<T>
		{
			result.SetTo(marginalB);
			return result;
		}

		/// <summary>
		/// VMP message to 'def'
		/// </summary>
		/// <param name="marginal">Current 'marginal'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'def' conditioned on the given values.
		/// </para></remarks>
		public static T DefAverageLogarithm<T>([IsReturned] T marginalB, T result)
			where T : SettableTo<T>
		{
			result.SetTo(marginalB);
			return result;
		}
	}
	/// <summary>
	/// Provides outgoing VMP messages for <see cref="Factor.Variable&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Variable<>")]
	[Buffers("marginalB")]
	[Quality(QualityBand.Preview)]
	public static class VariableNoInitVmpBufferOp
	{
		[Skip]
		public static T MarginalBInit<T>([IgnoreDependency] T def)
			where T : ICloneable
		{
			return (T)def.Clone();
		}
	}

	/// <summary>
	/// Provides outgoing VMP messages for <see cref="Factor.Variable&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "VariableInit<>", Default=true)]
	[Buffers("marginalB")]
	[Quality(QualityBand.Preview)]
	public static class VariableInitVmpBufferOp
	{
		// def is included to get its type as a constraint.  not needed if we could bind on return type.
		public static T MarginalBInit<T>([IgnoreDependency] T def, [SkipIfUniform] T init)
			where T : ICloneable
		{
			return (T)init.Clone();
		}
	}
#endif

	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.DerivedVariable{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "DerivedVariable<>")]
	[FactorMethod(typeof(Factor), "DerivedVariableInit<>")]
	[Quality(QualityBand.Preview)]
	public static class DerivedVariableOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(use,def,marginal))</c>.
		/// </para></remarks>
		[Skip]
		public static double LogAverageFactor()
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(use,def,marginal))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio()
		{
			return 0.0;
		}

		/// <summary>
		/// EP message to 'marginal'
		/// </summary>
		/// <param name="Use">Incoming message from 'use'.</param>
		/// <param name="Def">Incoming message from 'def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'marginal' as the random arguments are varied.
		/// The formula is <c>proj[p(marginal) sum_(use,def) p(use,def) factor(use,def,marginal)]/p(marginal)</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		[MultiplyAll]
		public static T MarginalAverageConditional<T>(T Use, T Def, T result)
				where T : SettableToProduct<T>
		{
			result.SetToProduct(Def, Use);
			return result;
		}
		[Skip]
		public static T MarginalAverageConditionalInit<T>([IgnoreDependency] T def)
			where T : ICloneable
		{
			return (T)def.Clone();
		}

		/// <summary>
		/// EP message to 'use'
		/// </summary>
		/// <param name="Def">Incoming message from 'def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'use' as the random arguments are varied.
		/// The formula is <c>proj[p(use) sum_(def) p(def) factor(use,def,marginal)]/p(use)</c>.
		/// </para></remarks>
		public static T UseAverageConditional<T>([IsReturned] T Def)
		{
			return Def;
		}

		/// <summary>
		/// EP message to 'def'
		/// </summary>
		/// <param name="Use">Incoming message from 'use'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'def' as the random arguments are varied.
		/// The formula is <c>proj[p(def) sum_(use) p(use) factor(use,def,marginal)]/p(def)</c>.
		/// </para></remarks>
		public static T DefAverageConditional<T>([IsReturned] T Use)
		{
			return Use;
		}
	}
	[FactorMethod(typeof(Factor), "DerivedVariableGibbs<>")]
	[FactorMethod(typeof(Factor), "DerivedVariableInitGibbs<>")]
	[Quality(QualityBand.Preview)]
	public static class DerivedVariableGibbsOp
	{
		#region Gibbs messages
		/// <summary>
		/// Evidence message for Gibbs.
		/// </summary>
		[Skip]
		public static double GibbsEvidence()
		{
			return 0.0;
		}

		/// <summary>
		/// Gibbs sample message to 'Uses'
		/// </summary>
		/// <typeparam name="TMarginalDist">Gibbs marginal type</typeparam>
		/// <typeparam name="TDef">Definition distribution type</typeparam>
		/// <typeparam name="T">Domain type</typeparam>
		/// <param name="marginal">The Gibbs marginal</param>
		/// <param name="resultIndex">'Uses' index for result (unused)</param>
		/// <param name="result">Result</param>
		/// <returns></returns>
		/// <remarks><para>
		/// The outgoing message is the current Gibbs sample.
		/// </para></remarks>
		public static T UseGibbs<TDist, T>([SkipIfUniform] GibbsMarginal<TDist, T> marginal, TDist def, T result)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			// This method must depend on Def, even though Def isn't used, in order to get the right triggers
			return marginal.LastSample;
		}

		/// <summary>
		/// Gibbs sample message to 'Uses'
		/// </summary>
		/// <typeparam name="T">Domain type</typeparam>
		/// <param name="def"></param>
		/// <param name="resultIndex">'Uses' index for result (unused)</param>
		/// <param name="result">Result</param>
		/// <returns></returns>
		/// <remarks><para>
		/// The outgoing message is the current Gibbs sample.
		/// </para></remarks>
		public static T UseGibbs<T>([IsReturned] T def, T result)
		{
			return def;
		}

		/// <summary>
		/// Gibbs distribution message to 'Def'
		/// </summary>
		/// <typeparam name="TDist">Distribution type</typeparam>
		/// <typeparam name="T">Domain type</typeparam>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the product of all the 'Uses' messages.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		public static TDist DefGibbs<TDist, T>([IsReturned] TDist Use)
		{
			return Use;
		}
		public static T DefGibbs<TDist, T>([SkipIfUniform] GibbsMarginal<TDist, T> marginal, T result)
			 where TDist : IDistribution<T>, Sampleable<T>
		{
			return marginal.LastSample;
		}

		/// <summary>
		/// Gibbs message to 'Marginal' for distribution Def
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="to_marginal">Previous outgoing message to 'marginal'.</param>
		/// <returns><paramref name="to_marginal"/></returns>
		/// <remarks><para>
		/// The outgoing message is the product of 'Def' and 'Uses' messages.
		/// </para></remarks>
		[Stochastic]
		[SkipIfAllUniform]
		public static GibbsMarginal<TDist, T> MarginalGibbs<TDist, T>(
			TDist Use,
			[SkipIfUniform]TDist Def,
			GibbsMarginal<TDist, T> to_marginal)
			where TDist : IDistribution<T>, SettableToProduct<TDist>, Sampleable<T>
		{
			GibbsMarginal<TDist, T> result = to_marginal;
			TDist marginal = result.LastConditional;
			marginal.SetToProduct(Def, Use);
			result.LastConditional = marginal;
			// Allow a sample to be drawn from the last conditional, and add it to the sample
			// list and conditional list
			result.PostUpdate();
			return result;
		}

		/// <summary>
		/// Gibbs message to 'Marginal' for sample Def
		/// </summary>
		/// <typeparam name="TDist"></typeparam>
		/// <param name="Def"></param>
		/// <param name="to_marginal">Previous outgoing message to 'marginal'.</param>
		/// <returns><paramref name="to_marginal"/></returns>
		[Stochastic] // must be labelled Stochastic to get correct schedule, even though it isn't Stochastic
		public static GibbsMarginal<TDist, T> MarginalGibbs<TDist, T>(
			T Def,
			GibbsMarginal<TDist, T> to_marginal)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			GibbsMarginal<TDist, T> result = to_marginal;
			TDist marginal = result.LastConditional;
			marginal.Point = Def;
			result.LastConditional = marginal;
			// Allow a sample to be drawn from the last conditional, and add it to the sample
			// list and conditional list
			result.PostUpdate();
			return result;
		}
		/// <summary>
		/// Gibbs message to 'Marginal' for sample Use
		/// </summary>
		/// <typeparam name="TDist"></typeparam>
		/// <param name="Use"></param>
		/// <param name="to_marginal">Previous outgoing message to 'marginal'.</param>
		/// <returns><paramref name="to_marginal"/></returns>
		[Stochastic] // must be labelled Stochastic to get correct schedule, even though it isn't Stochastic
		public static GibbsMarginal<TDist, T> MarginalGibbs<TDist, T>(
			T Use, [IgnoreDependency] TDist Def,
			GibbsMarginal<TDist, T> to_marginal)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			GibbsMarginal<TDist, T> result = to_marginal;
			TDist marginal = result.LastConditional;
			marginal.Point = Use;
			result.LastConditional = marginal;
			// Allow a sample to be drawn from the last conditional, and add it to the sample
			// list and conditional list
			result.PostUpdate();
			return result;
		}

		#endregion
	}
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.DerivedVariable{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "DerivedVariableVmp<>")]
	[FactorMethod(typeof(Factor), "DerivedVariableInitVmp<>")]
	[Quality(QualityBand.Preview)]
	public static class DerivedVariableVmpOp
	{
		// Derived variables send no evidence messages.
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(use,def,marginal))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'marginal'
		/// </summary>
		/// <param name="Def">Incoming message from 'def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'marginal'.
		/// The formula is <c>exp(sum_(def) p(def) log(factor(use,def,marginal)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		public static T MarginalAverageLogarithm<T, TDef>([IsReturned] TDef Def, T result)
		 where T : SettableTo<TDef>
		{
			result.SetTo(Def);
			return result;
		}
		[Skip]
		public static T MarginalAverageLogarithmInit<T>([IgnoreDependency] T def)
			where T : ICloneable
		{
			return (T)def.Clone();
		}

		/// <summary>
		/// VMP message to 'use'
		/// </summary>
		/// <param name="Def">Incoming message from 'def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'use' as the random arguments are varied.
		/// The formula is <c>proj[sum_(def) p(def) factor(use,def,marginal)]</c>.
		/// </para></remarks>
		public static T UseAverageLogarithm<T, TDef>([IsReturned] TDef Def, T result)
		 where T : SettableTo<TDef>
		{
			result.SetTo(Def);
			return result;
		}
		[Skip]
		public static T UseAverageLogarithmInit<T>([IgnoreDependency] T def)
			where T : ICloneable
		{
			return (T)def.Clone();
		}

		// must have upward Trigger to match the Trigger on UsesEqualDef.UsesAverageLogarithm
		/// <summary>
		/// VMP message to 'def'
		/// </summary>
		/// <param name="Use">Incoming message from 'use'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'def' with 'use' integrated out.
		/// The formula is <c>sum_use p(use) factor(use,def,marginal)</c>.
		/// </para></remarks>
		public static T DefAverageLogarithm<T>([IsReturned] T Use, T result)
		where T : SettableTo<T>
		{
			result.SetTo(Use);
			return result;
		}
	}
}
