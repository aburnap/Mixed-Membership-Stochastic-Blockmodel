// (C) Copyright 2008 Microsoft Research Cambridge
#define SpecializeArrays
#define MinimalGenericTypeParameters

using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;
using MicrosoftResearch.Infer.Collections;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Messages from random variables that are defined by deterministic factors.
	/// </summary>
	[FactorMethod(typeof(Factor), "Replicate<>", Default = true)]
	[Buffers("marginal", "toDef")]
	[Quality(QualityBand.Preview)]
	public static class ReplicateOp_Divide
	{
		public static T DefAverageConditional<T>([IsReturned] T toDef, T result)
			where T : SettableTo<T>
		{
			result.SetTo(toDef);
			return result;
		}

		/// <summary>
		/// EP message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="marginal">Buffer 'marginal'.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Uses' conditioned on the given values.
		/// </para></remarks>
		//[SkipIfAllUniform]
		// Uses dependency must be ignored for Sequential schedule
		// Uses is marked ignore because the forward message does not really depend on the backward message
		public static T UsesAverageConditional<T>([Indexed, Cancels] T Uses, [Fresh, SkipIfUniform] T marginal, int resultIndex, T result)
			where T : SettableToRatio<T>
		{
			result.SetToRatio(marginal, Uses);
			return result;
		}

		/// <summary>
		/// Initialise the buffer 'marginal'
		/// </summary>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <returns>Initial value of buffer 'marginal'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip] // this is needed to instruct the scheduler to treat marginal as uninitialized
		public static T MarginalInit<T>([SkipIfUniform] T Def)
				where T : ICloneable, SettableToUniform
		{
			return ArrayHelper.MakeUniform(Def);
		}

		/// <summary>
		/// Update the buffer 'marginal'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[SkipIfAllUniform]
		[MultiplyAll]
		public static T Marginal<T>([Fresh] T toDef, T Def, T result)
				where T : SettableToProduct<T>
		{
			result.SetToProduct(Def, toDef);
			return result;
		}

		public static T MarginalIncrement<T>(T result, [SkipIfUniform] T use, [SkipIfUniform] T def)
			where T : SettableToProduct<T>
		{
			result.SetToProduct(use, def);
			return result;
		}

		/// <summary>
		/// Initialise the buffer 'toDef'
		/// </summary>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <returns>Initial value of buffer 'toDef'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip] // this is needed to instruct the scheduler to treat the buffer as uninitialized
		public static T ToDefInit<T>(T Def)
			where T : ICloneable, SettableToUniform
		{
			// must construct from Def instead of Uses because Uses array may be empty
			return ArrayHelper.MakeUniform(Def);
		}

		[SkipIfAllUniform]
		[MultiplyAll]
		public static T ToDef<T>(IList<T> Uses, T result)
			where T : SettableToProduct<T>, SettableTo<T>, SettableToUniform
		{
			return Distribution.SetToProductOfAll(result, Uses);
		}
	}

	/// <summary>
	/// Messages from random variables that are defined by deterministic factors.
	/// </summary>
	[FactorMethod(typeof(Factor), "Replicate<>", Default = false)]
	[Buffers("marginal")]
	[Quality(QualityBand.Preview)]
	public static class Replicate2BufferOp
	{
#if false
		/// <summary>
		/// EP message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="marginal">Buffer 'marginal'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Uses' conditioned on the given values.
		/// </para></remarks>
		//[SkipIfAllUniform]
		public static TListRet UsesAverageConditional<T, TList, TListRet>(TList Uses, [Fresh, SkipIfUniform] T marginal, TListRet result)
			where TList : IList<T>
			where TListRet : IList<T>
			where T : SettableToRatio<T>
		{
			for (int i = 0; i < result.Count; i++) {
				T dist = result[i];
				dist.SetToRatio(marginal, Uses[i]);
				result[i] = dist;
			}
			return result;
		}
#endif

		/// <summary>
		/// EP message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="marginal">Buffer 'marginal'.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Uses' conditioned on the given values.
		/// </para></remarks>
		//[SkipIfAllUniform]
		// Uses dependency must be ignored for Sequential schedule
		public static T UsesAverageConditional<T>([MatchingIndex, IgnoreDependency] IList<T> Uses, [IgnoreDependency, SkipIfUniform] T Def, [Fresh, SkipIfUniform] T marginal, int resultIndex, T result)
			where T : SettableToRatio<T>, SettableToProduct<T>, SettableTo<T>
		{
			if (resultIndex < 0 || resultIndex >= Uses.Count) throw new ArgumentOutOfRangeException("resultIndex");
			if (Uses.Count == 1) { result.SetTo(Def); return result; }
			if (true) {
				try {
					result.SetToRatio(marginal, Uses[resultIndex]);
					return result;
				} catch (DivideByZeroException) {
					return ReplicateOp_NoDivide.UsesAverageConditional(Uses, Def, resultIndex, result);
				}
			} else {
				// check that ratio is same as product
				result.SetToRatio(marginal, Uses[resultIndex]);
				T result2 = (T)((ICloneable)result).Clone();
				ReplicateOp_NoDivide.UsesAverageConditional(Uses, Def, resultIndex, result2);
				double err = ((Diffable)result).MaxDiff(result2);
				if (err > 1e-4) Console.WriteLine(err);
				return result;
			}
		}
#if SpecializeArrays
		/// <summary>
		/// EP message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="marginal">Buffer 'marginal'.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Uses' conditioned on the given values.
		/// </para></remarks>
		//[SkipIfAllUniform]
		// Uses dependency must be ignored for Sequential schedule
		public static T UsesAverageConditional<T>([MatchingIndex, IgnoreDependency] T[] Uses, [IgnoreDependency, SkipIfUniform] T Def, [Fresh, SkipIfUniform] T marginal, int resultIndex, T result)
			where T : SettableToRatio<T>, SettableToProduct<T>, SettableTo<T>
		{
			if (resultIndex < 0 || resultIndex >= Uses.Length) throw new ArgumentOutOfRangeException("resultIndex");
			if (Uses.Length == 1) { result.SetTo(Def); return result; }
			try {
				result.SetToRatio(marginal, Uses[resultIndex]);
				return result;
			} catch (DivideByZeroException) {
				return ReplicateOp_NoDivide.UsesAverageConditional(Uses, Def, resultIndex, result);
			}
		}
#endif
		[Skip]
		public static ArrayType UsesAverageConditionalInit<T, ArrayType>([IgnoreDependency] T Def, int count, IArrayFactory<T, ArrayType> factory)
			 where T : ICloneable
		{
			return factory.CreateArray(count, i => (T)Def.Clone());
		}

		/// <summary>
		/// Initialise the buffer 'marginal'
		/// </summary>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <returns>Initial value of buffer 'marginal'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip] // this is needed to instruct the scheduler to treat marginal as uninitialized
		public static T MarginalInit<T>([SkipIfUniform] T Def)
				where T : ICloneable
		{
			return (T)Def.Clone();
		}

		/// <summary>
		/// Update the buffer 'marginal'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[SkipIfAllUniform]
		[MultiplyAll]
		public static T Marginal<T>(IList<T> Uses, T Def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			return ReplicateOp_NoDivide.MarginalAverageConditional(Uses, Def, result);
		}

		public static T MarginalIncrement<T>(T result, [SkipIfUniform] T use, [SkipIfUniform] T def)
			where T : SettableToProduct<T>
		{
			result.SetToProduct(use, def);
			return result;
		}

#if MinimalGenericTypeParameters
		/// <summary>
		/// EP message to 'Def'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Def' as the random arguments are varied.
		/// The formula is <c>proj[p(Def) sum_(Uses) p(Uses) factor(Uses,Def,Marginal)]/p(Def)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		[MultiplyAll]
		public static T DefAverageConditional<T>([SkipIfAllUniform] IList<T> Uses, T result)
			where T : SettableToProduct<T>, SettableTo<T>, SettableToUniform
		{
			return Distribution.SetToProductOfAll(result, Uses);
		}
#if SpecializeArrays
		/// <summary>
		/// EP message to 'Def'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Def' as the random arguments are varied.
		/// The formula is <c>proj[p(Def) sum_(Uses) p(Uses) factor(Uses,Def,Marginal)]/p(Def)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		[MultiplyAll]
		public static T DefAverageConditional<T>([SkipIfAllUniform] T[] Uses, T result)
		where T : SettableToProduct<T>, SettableTo<T>, SettableToUniform
		{
			return Distribution.SetToProductOfAll(result, Uses);
		}
#endif
#else
			public static T DefAverageConditional<T,TUses>([SkipIfAllUniform] IList<TUses> Uses, T result)
            where T : SettableToProduct<TUses>, SettableTo<TUses>, TUses, SettableToUniform
        {
            return Distribution.SetToProductOfAll(result, Uses);
        }
#if SpecializeArrays
        public static T DefAverageConditional<T,TUses>([SkipIfAllUniform] TUses[] Uses, T result)
            where T : SettableToProduct<TUses>, SettableTo<TUses>, TUses, SettableToUniform
        {
            return Distribution.SetToProductOfAll(result, Uses);
        }
#endif
#endif
	}
	/// <summary>
	/// Messages from random variables that are defined by deterministic factors.
	/// </summary>
	[FactorMethod(typeof(Factor), "ReplicateWithMarginal<>", Default = true)]
	[Quality(QualityBand.Preview)]
	public static class ReplicateBufferOp
	{
		/// <summary>
		/// EP message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="to_marginal">Previous outgoing message to 'marginal'.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Uses' conditioned on the given values.
		/// </para></remarks>
		//[SkipIfAllUniform]
		public static T UsesAverageConditional<T>([AllExceptIndex] IList<T> Uses, [SkipIfUniform] T Def, [SkipIfUniform, Fresh] T to_marginal, int resultIndex, T result)
			where T : SettableToRatio<T>, SettableToProduct<T>, SettableTo<T>
		{
			if (resultIndex < 0 || resultIndex >= Uses.Count) throw new ArgumentOutOfRangeException("resultIndex");
			if (Uses.Count == 1) { result.SetTo(Def); return result; }
			try {
				result.SetToRatio(to_marginal, Uses[resultIndex]);
				return result;
			} catch (DivideByZeroException) {
				return ReplicateOp_NoDivide.UsesAverageConditional(Uses, Def, resultIndex, result);
			}
		}
#if SpecializeArrays
		/// <summary>
		/// EP message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="to_marginal">Previous outgoing message to 'marginal'.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Uses' conditioned on the given values.
		/// </para></remarks>
		//[SkipIfAllUniform]
		public static T UsesAverageConditional<T>([AllExceptIndex] T[] Uses, [SkipIfUniform] T Def, [SkipIfUniform, Fresh] T to_marginal, int resultIndex, T result)
			where T : SettableToRatio<T>, SettableToProduct<T>, SettableTo<T>
		{
			if (resultIndex < 0 || resultIndex >= Uses.Length) throw new ArgumentOutOfRangeException("resultIndex");
			if (Uses.Length == 1) { result.SetTo(Def); return result; }
			try {
				result.SetToRatio(to_marginal, Uses[resultIndex]);
				return result;
			} catch (DivideByZeroException) {
				return ReplicateOp_NoDivide.UsesAverageConditional(Uses, Def, resultIndex, result);
			}
		}
#endif
	}

	/// <summary>
	/// Messages from random variables that are defined by deterministic factors.
	/// </summary>
	[FactorMethod(typeof(Factor), "Replicate<>", Default = true)]
	[FactorMethod(typeof(Factor), "ReplicateWithMarginal<>", Default = true)]
	[Quality(QualityBand.Mature)]
	public static class ReplicateGibbsOp<T>
	{
		/// <summary>
		/// Gibbs evidence
		/// </summary>
		/// <returns></returns>
		/// <remarks><para>
		/// Returns 0.0
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio<TDist>([SkipIfAllUniform] IList<TDist> Uses, T Def)
			 where TDist : IDistribution<T>
		{
			return 0.0;
		}

		/// <summary>
		/// Gibbs evidence
		/// </summary>
		/// <returns></returns>
		/// <remarks><para>
		/// Returns 0.0
		/// </para></remarks>
		[Skip]
		public static double GibbsEvidence<TDist>(IList<TDist> Uses, T Def)
			where TDist : IDistribution<T>
		{
			return 0.0;
		}

#if false
		public static T UsesGibbs<TDist>([SkipIfUniform] GibbsMarginal<TDist, T> marginal, int resultIndex, T result)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			return marginal.LastSample;
		}
#elif false
		public static T UsesGibbs<TDist>([SkipIfUniform] GibbsMarginal<TDist, T> marginal, TDist def, int resultIndex, T result)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			return marginal.LastSample;
		}
		public static T UsesGibbs<TDist>([SkipIfUniform] GibbsMarginal<TDist, T> marginal, T def, int resultIndex, T result)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			if (def is bool[]) {
				if (!Util.ValueEquals((bool[])(object)def, (bool[])(object)marginal.LastSample)) throw new Exception("gotcha");
			}
			else if (def is double[]) {
				if (!Util.ValueEquals((double[])(object)def, (double[])(object)marginal.LastSample)) throw new Exception("gotcha");
			} else if (!def.Equals(marginal.LastSample)) throw new Exception("gotcha");
			return marginal.LastSample;
		}
#else
		/// <summary>
		/// Gibbs sample message to 'Uses'
		/// </summary>
		/// <param name="to_marginal">Incoming message from 'Marginal'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the current Gibbs sample.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="to_marginal"/> is not a proper distribution</exception>
		public static T UsesGibbs<TDist>([SkipIfUniform] GibbsMarginal<TDist, T> to_marginal, TDist def, int resultIndex, T result)
			 where TDist : IDistribution<T>, Sampleable<T>
		{
			// This method must depend on Def, even though Def isn't used, in order to get the right triggers
			return to_marginal.LastSample;
		}

		/// <summary>
		/// Gibbs sample message to 'Uses'
		/// </summary>
		/// <param name="def">Incoming message from 'Def'.</param>
		/// <param name="resultIndex">'Uses' index for result (unused)</param>
		/// <param name="result">Result</param>
		/// <returns></returns>
		/// <remarks><para>
		/// The outgoing message is the message from def.
		/// </para></remarks>
		public static T UsesGibbs([IsReturned] T def, int resultIndex, T result)
		{
			return def;
		}
#endif

#if true // until .NET 4
		/// <summary>
		/// Gibbs distribution message to 'Uses'
		/// </summary>
		/// <typeparam name="TDist">Distribution type</typeparam>
		/// <param name="def">Incoming message from 'Def'.</param>
		/// <param name="resultIndex">'Uses' index for result (unused)</param>
		/// <param name="result">Result</param>
		/// <returns></returns>
		/// <remarks><para>
		/// The outgoing message is the message from def.
		/// </para></remarks>
		public static TDist UsesGibbs<TDist>([IsReturned] TDist def, int resultIndex, TDist result)
			where TDist : IDistribution<T>
		{
			return def;
		}
#endif

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
		[MultiplyAll]
		public static TDist DefGibbs<TDist>(
			[SkipIfAllUniform]IList<TDist> Uses,
			TDist result)
			where TDist : IDistribution<T>, Sampleable<T>, SettableTo<TDist>, SettableToProduct<TDist>, SettableToRatio<TDist>
		{
			return ReplicateOp_NoDivide.DefAverageConditional(Uses, result);
		}
#if SpecializeArrays
		/// <summary>
		/// Gibbs distribution message to 'Def'
		/// </summary>
		/// <typeparam name="TDist">Distribution type</typeparam>
		/// <typeparam name="T">Domain type</typeparam>
		/// <param name="Uses">Uses</param>
		/// <param name="result">Result</param>
		/// <returns></returns>
		[MultiplyAll]
		public static TDist DefGibbs<TDist>(
			[SkipIfAllUniform]TDist[] Uses,
			TDist result)
			where TDist : IDistribution<T>, Sampleable<T>, SettableToProduct<TDist>, SettableTo<TDist>
		{
			return ReplicateOp_NoDivide.DefAverageConditional(Uses, result);
		}
#endif

		/// <summary>
		/// Gibbs sample message to 'Def'
		/// </summary>
		/// <typeparam name="TDist">Distribution type</typeparam>
		/// <typeparam name="T">Domain type</typeparam>
		/// <param name="to_marginal">The Gibbs marginal</param>
		/// <param name="result">Result</param>
		/// <returns></returns>
		/// <remarks><para>
		/// The outgoing message is the current Gibbs sample.
		/// </para></remarks>
		public static T DefGibbs<TDist>([SkipIfUniform] GibbsMarginal<TDist, T> to_marginal, T result)
			 where TDist : IDistribution<T>, Sampleable<T>
		{
			return to_marginal.LastSample;
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
		public static GibbsMarginal<TDist, T> MarginalGibbs<TDist>(
			IList<TDist> Uses,
			[SkipIfUniform]TDist Def,
			GibbsMarginal<TDist, T> to_marginal)
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToRatio<TDist>, SettableTo<TDist>, Sampleable<T>
		{
			GibbsMarginal<TDist, T> result = to_marginal;
			TDist marginal = result.LastConditional;
			marginal.SetTo(Def);
			marginal = Distribution.SetToProductWithAll(marginal, Uses);
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
		public static GibbsMarginal<TDist, T> MarginalGibbs<TDist>(
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
		/// Gibbs message to 'Marginal' for sample Uses
		/// </summary>
		/// <typeparam name="TDist"></typeparam>
		/// <param name="Uses"></param>
		/// <param name="to_marginal">Previous outgoing message to 'marginal'.</param>
		/// <returns><paramref name="to_marginal"/></returns>
		[Stochastic] // must be labelled Stochastic to get correct schedule, even though it isn't Stochastic
		public static GibbsMarginal<TDist, T> MarginalGibbs<TDist>(
			T[] Uses,
			GibbsMarginal<TDist, T> to_marginal)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			GibbsMarginal<TDist, T> result = to_marginal;
			TDist marginal = result.LastConditional;
			if (Uses.Length != 1) throw new ArgumentException("Uses.Length ("+Uses.Length+") != 1");
			marginal.Point = Uses[0];
			result.LastConditional = marginal;
			// Allow a sample to be drawn from the last conditional, and add it to the sample
			// list and conditional list
			result.PostUpdate();
			return result;
		}
		[Skip]
		public static GibbsMarginal<TDist, T> MarginalGibbsInit<TDist>([IgnoreDependency] TDist def)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			return new GibbsMarginal<TDist, T>(def, 100, 1, true, true, true);
		}
	}

	/// <summary>
	/// Provides outgoing Gibbs messages for <see cref="Factor.UsesEqualDef&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "ReplicateWithMarginalGibbs<>")]
	[Buffers("sample", "conditional", "marginalEstimator", "sampleAcc", "conditionalAcc")]
	[Quality(QualityBand.Mature)]
	public static class ReplicateGibbsOp2<T>
	{
		[Skip]
		public static TDist ConditionalInit<TDist>([IgnoreDependency] TDist to_marginal)
			where TDist : ICloneable
		{
			return (TDist)to_marginal.Clone();
		}
		public static TDist Conditional<TDist>(T Def, TDist result)
			where TDist : HasPoint<T>
		{
			result.Point = Def;
			return result;
		}
		public static TDist Conditional<TDist>(IList<TDist> Uses, [SkipIfAnyUniform]TDist Def, TDist result)
			where TDist : SettableTo<TDist>, SettableToProduct<TDist>
		{
			result.SetTo(Def);
			result = Distribution.SetToProductWithAll(result, Uses);
			return result;
		}
		[Stochastic]
		public static T Sample<TDist>([IgnoreDependency] TDist to_marginal, [Proper]TDist conditional)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			return conditional.Sample();
		}

		public static BurnInAccumulator<TDist> MarginalEstimatorInit<TDist>([IgnoreDependency] TDist to_marginal, int burnIn)
			 where TDist : IDistribution<T>
		{
			Accumulator<TDist> est = (Accumulator<TDist>)ArrayEstimator.CreateEstimator<TDist, T>(to_marginal, true);
			return new BurnInAccumulator<TDist>(burnIn, 1, est);
		}
		public static TAcc MarginalEstimator<TDist, TAcc>([Proper]TDist conditional, TAcc marginalEstimator)
			 where TAcc : Accumulator<TDist>
		{
			marginalEstimator.Add(conditional);
			return marginalEstimator;
		}

		public static TDist MarginalGibbs<TDist>(BurnInAccumulator<TDist> marginalEstimator, TDist result)
		{
			return ((Estimator<TDist>)marginalEstimator.Accumulator).GetDistribution(result);
		}

		public static Accumulator<T> SampleAccInit(ICollection<T> to_samples, int burnIn, int thin)
		{
			return new BurnInAccumulator<T>(burnIn, thin, new AccumulateIntoCollection<T>(to_samples));
		}
		public static Accumulator<T> SampleAcc(T sample, Accumulator<T> sampleAcc)
		{
			sampleAcc.Add(sample);
			return sampleAcc;
		}
		public static TList SamplesGibbs<TList>(Accumulator<T> sampleAcc, TList result)
			where TList : ICollection<T>
		{
			// do nothing since result was already modified by sampleAcc
			return result;
		}

		public static Accumulator<TDist> ConditionalAccInit<TDist>(ICollection<TDist> to_conditionals, int burnIn, int thin)
		{
			return new BurnInAccumulator<TDist>(burnIn, thin, new AccumulateIntoCollection<TDist>(to_conditionals));
		}
		public static Accumulator<TDist> ConditionalAcc<TDist>(TDist conditional, Accumulator<TDist> conditionalAcc)
			where TDist : ICloneable
		{
			conditionalAcc.Add((TDist)conditional.Clone());
			return conditionalAcc;
		}
		public static TDistList ConditionalsGibbs<TDist, TDistList>(Accumulator<TDist> conditionalAcc, TDistList result)
			where TDistList : ICollection<TDist>
		{
			// do nothing since result was already modified by Acc
			return result;
		}

		/// <summary>
		/// Gibbs evidence
		/// </summary>
		/// <returns></returns>
		/// <remarks><para>
		/// Returns 0.0
		/// </para></remarks>
		[Skip]
		public static double GibbsEvidence<TDist>(IList<TDist> Uses, T Def)
			where TDist : IDistribution<T>
		{
			return 0.0;
		}

		/// <summary>
		/// Gibbs sample message to 'Uses'
		/// </summary>
		/// <param name="def">Incoming message from 'Def'.</param>
		/// <param name="resultIndex">'Uses' index for result (unused)</param>
		/// <param name="result">Result</param>
		/// <returns></returns>
		/// <remarks><para>
		/// The outgoing message is the message from def.
		/// </para></remarks>
		public static T UsesGibbs([IsReturned] T def, int resultIndex, T result)
		{
			return def;
		}
		/// <summary>
		/// Gibbs sample message to 'Uses'
		/// </summary>
		/// <param name="sample">Current sample.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the current Gibbs sample.
		/// </para></remarks>
		public static T UsesGibbs<TDist>(TDist def, T sample, int resultIndex, T result)
			where TDist : IDistribution<T>
		{
			// This method must depend on Def, even though Def isn't used, in order to get the right triggers
			return sample;
		}

		/// <summary>
		/// Gibbs distribution message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the product of the 'Def' message with all 'Uses' messages except the current
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		public static TDist UsesGibbs<TDist>(
			[IgnoreDependency] ICollection<TDist> Uses,
			[IsReturned] TDist Def,
			int resultIndex, TDist result)
			where TDist : IDistribution<T>, Sampleable<T>, SettableTo<TDist>, SettableToProduct<TDist>, SettableToRatio<TDist>
		{
			if (resultIndex < 0 || resultIndex >= Uses.Count) throw new ArgumentOutOfRangeException("resultIndex");
			if (Uses.Count > 1) throw new ArgumentException("Uses.Count > 1");
			result.SetTo(Def);
			return result;
		}

		[Skip]
		public static T UsesGibbsInit<T>([IgnoreDependency] T Def, int resultIndex)
			where T : ICloneable
		{
			return (T)Def.Clone();
		}

		public static T DefGibbs<TDist>(TDist def, [IsReturned] T sample)
			 where TDist : IDistribution<T>
		{
			return sample;
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
		//[MultiplyAll]
		public static TDist DefGibbs<TDist>(
			[SkipIfAllUniform]IList<TDist> Uses,
			TDist result)
			where TDist : IDistribution<T>, Sampleable<T>, SettableTo<TDist>, SettableToProduct<TDist>, SettableToRatio<TDist>
		{
			result.SetToUniform();
			result = Distribution.SetToProductWithAll(result, Uses);
			return result;
		}
	}

	/// <summary>
	/// Messages from random variables that are defined by deterministic factors.
	/// </summary>
	[FactorMethod(typeof(Factor), "Replicate<>", Default = false)]
	[FactorMethod(typeof(Factor), "ReplicateWithMarginal<>", Default = false)]
	[Quality(QualityBand.Mature)]
	public static class ReplicateOp_NoDivide
	{
		/// <summary>
		/// EP message to 'Marginal'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Marginal' as the random arguments are varied.
		/// The formula is <c>proj[p(Marginal) sum_(Uses,Def) p(Uses,Def) factor(Uses,Def,Marginal)]/p(Marginal)</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		[MultiplyAll]
		public static T MarginalAverageConditional<T>(IList<T> Uses, T Def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			result.SetTo(Def);
			return Distribution.SetToProductWithAll(result, Uses);
		}
#if SpecializeArrays
		/// <summary>
		/// EP message to 'Marginal'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Marginal' as the random arguments are varied.
		/// The formula is <c>proj[p(Marginal) sum_(Uses,Def) p(Uses,Def) factor(Uses,Def,Marginal)]/p(Marginal)</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		[MultiplyAll]
		public static T MarginalAverageConditional<T>(T[] Uses, T Def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			result.SetTo(Def);
			return Distribution.SetToProductWithAll(result, Uses);
		}
#endif

#if false
		/// <summary>
		/// EP message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Uses' as the random arguments are varied.
		/// The formula is <c>proj[p(Uses) sum_(Def) p(Def) factor(Uses,Def,Marginal)]/p(Uses)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		//[SkipIfAllUniform]
		public static TListRet UsesAverageConditional<T, TList, TListRet>(TList Uses, [SkipIfUniform] T Def, TListRet result)
			where TList : IList<T>
			where TListRet : IList<T>
			where T : SettableToProduct<T>, SettableTo<T>, SettableToRatio<T>
		{
			if (Uses.Count < 5) {
				for (int i = 0; i < result.Count; i++) {
					T dist = result[i];
					dist.SetTo(Def);
					dist.SetTo(Distribution.SetToProductWithAllExcept(dist, Uses, i));
					result[i] = dist;
				}
			} else if (true) {
				T productOfAll = (T)((ICloneable)Def).Clone();
				productOfAll = Distribution.SetToProductWithAll(productOfAll, Uses);
				for (int i = 0; i < result.Count; i++) {
					T dist = result[i];
					dist.SetToRatio(productOfAll, Uses[i]);
					result[i] = dist;
				}
			} else { // for debugging
				for (int i = 0; i < result.Count; i++) {
					T dist = result[i];
					dist.SetTo(Def);
					dist.SetTo(Distribution.SetToProductWithAllExcept(dist, Uses, i));
					result[i] = dist;
				}

				T productOfAll = (T)((ICloneable)Def).Clone();
				productOfAll = Distribution.SetToProductWithAll(productOfAll, Uses);
				for (int i = 0; i < result.Count; i++) {
					T dist = result[i];
					T olddist = (T)((ICloneable)dist).Clone();
					dist.SetTo(Def);
					dist.SetToRatio(productOfAll, Uses[i]);
					double maxdiff = ((Diffable)olddist).MaxDiff(dist);
					if (maxdiff > 1e-10) throw new OverflowException();
					// for SparseGPs where Uses[i] is a rank1Pot not on a basis point, division gives different answer than multiplication
					//Console.WriteLine("[{1}] maxdiff = {0}", maxdiff, StringUtil.TypeToString(typeof(T)));
					result[i] = dist;
				}
			}
			return result;
		}
#endif

		/// <summary>
		/// EP message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Uses' as the random arguments are varied.
		/// The formula is <c>proj[p(Uses) sum_(Def) p(Def) factor(Uses,Def,Marginal)]/p(Uses)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		//[SkipIfAllUniform]
		public static T UsesAverageConditional<T>([AllExceptIndex] IList<T> Uses, [SkipIfAllUniform] T Def, int resultIndex, T result)
			where T : SettableToProduct<T>, SettableTo<T>
		{
			if (resultIndex < 0 || resultIndex >= Uses.Count) throw new ArgumentOutOfRangeException("resultIndex");
			result.SetTo(Def);
			return Distribution.SetToProductWithAllExcept(result, Uses, resultIndex);
		}
#if SpecializeArrays
		/// <summary>
		/// EP message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Uses' as the random arguments are varied.
		/// The formula is <c>proj[p(Uses) sum_(Def) p(Def) factor(Uses,Def,Marginal)]/p(Uses)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		//[SkipIfAllUniform]
		public static T UsesAverageConditional<T>([AllExceptIndex] T[] Uses, [SkipIfAllUniform] T Def, int resultIndex, T result)
			where T : SettableToProduct<T>, SettableTo<T>
		{
			if (resultIndex < 0 || resultIndex >= Uses.Length) throw new ArgumentOutOfRangeException("resultIndex");
			result.SetTo(Def);
			return Distribution.SetToProductWithAllExcept(result, Uses, resultIndex);
		}
#endif

#if MinimalGenericTypeParameters
		/// <summary>
		/// EP message to 'Def'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Def' as the random arguments are varied.
		/// The formula is <c>proj[p(Def) sum_(Uses) p(Uses) factor(Uses,Def,Marginal)]/p(Def)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		[MultiplyAll]
		public static T DefAverageConditional<T>([SkipIfAllUniform] IList<T> Uses, T result)
			where T : SettableToProduct<T>, SettableTo<T>, SettableToUniform
		{
			return Distribution.SetToProductOfAll(result, Uses);
		}
#if SpecializeArrays
		/// <summary>
		/// EP message to 'Def'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Def' as the random arguments are varied.
		/// The formula is <c>proj[p(Def) sum_(Uses) p(Uses) factor(Uses,Def,Marginal)]/p(Def)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		[MultiplyAll]
		public static T DefAverageConditional<T>([SkipIfAllUniform] T[] Uses, T result)
		where T : SettableToProduct<T>, SettableTo<T>, SettableToUniform
		{
			return Distribution.SetToProductOfAll(result, Uses);
		}
#endif
#else
			public static T DefAverageConditional<T,TUses>([SkipIfAllUniform] IList<TUses> Uses, T result)
            where T : SettableToProduct<TUses>, SettableTo<TUses>, TUses, SettableToUniform
        {
            return Distribution.SetToProductOfAll(result, Uses);
        }
#if SpecializeArrays
        public static T DefAverageConditional<T,TUses>([SkipIfAllUniform] TUses[] Uses, T result)
            where T : SettableToProduct<TUses>, SettableTo<TUses>, TUses, SettableToUniform
        {
            return Distribution.SetToProductOfAll(result, Uses);
        }
#endif
#endif
	}

	/// <summary>
	/// Messages from random variables that are defined by deterministic factors.
	/// </summary>
	[FactorMethod(typeof(Factor), "Replicate<>", Default = true)]
	[FactorMethod(typeof(Factor), "ReplicateWithMarginal<>", Default = true)]
	[Quality(QualityBand.Mature)]
	public static class ReplicateOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Uses,Def,Marginal))</c>.
		/// </para></remarks>
		[Skip]
		public static double LogAverageFactor()
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Uses,Def) p(Uses,Def) factor(Uses,Def,Marginal) / sum_Uses p(Uses) messageTo(Uses))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio<T>([SkipIfAllUniform] IList<T> Uses, T Def, [Fresh] IList<T> to_Uses)
			where T : CanGetLogAverageOf<T>, SettableToProduct<T>, SettableTo<T>, ICloneable, SettableToUniform
		{
			return UsesEqualDefOp.LogEvidenceRatio(Uses, Def, to_Uses);
		}

		//-- VMP ----------------------------------------------------------------------------------------------

		// Deterministic variables send no evidence messages.
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Uses,Def,Marginal))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'Marginal'
		/// </summary>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'Marginal'.
		/// The formula is <c>exp(sum_(Def) p(Def) log(factor(Uses,Def,Marginal)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		public static T MarginalAverageLogarithm<T, TDef>([SkipIfAllUniform] TDef Def, T result)
		 where T : SettableTo<TDef>
		{
			return UsesAverageLogarithm<T, TDef>(Def, 0, result);
		}

		/// <summary>
		/// VMP message to 'Uses'
		/// </summary>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Uses' as the random arguments are varied.
		/// The formula is <c>proj[sum_(Def) p(Def) factor(Uses,Def,Marginal)]</c>.
		/// </para></remarks>
		public static T UsesAverageLogarithm<T, TDef>([IsReturned] TDef Def, int resultIndex, T result)
		 where T : SettableTo<TDef>
		{
			result.SetTo(Def);
			return result;
		}
		[Skip]
		public static T UsesDeriv<T>(T result)
			where T : SettableToUniform
		{
			result.SetToUniform();
			return result;
		}

		/// <summary>
		/// VMP message to 'Uses'
		/// </summary>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Uses' as the random arguments are varied.
		/// The formula is <c>proj[sum_(Def) p(Def) factor(Uses,Def,Marginal)]</c>.
		/// </para></remarks>
		public static T UsesAverageLogarithm<T, TDef>([IsReturnedInEveryElement] TDef Def, T result)
					where T : CanSetAllElementsTo<TDef>
		{
			result.SetAllElementsTo(Def);
			return result;
		}
		[Skip]
		public static ArrayType UsesInit<T, ArrayType>([IgnoreDependency] T Def, int count, IArrayFactory<T, ArrayType> factory)
			 where T : ICloneable
		{
			return factory.CreateArray(count, i => (T)Def.Clone());
		}

#if MinimalGenericTypeParameters
		/// <summary>
		/// VMP message to 'Def'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Def' with 'Uses' integrated out.
		/// The formula is <c>sum_Uses p(Uses) factor(Uses,Def,Marginal)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		[MultiplyAll]
		public static T DefAverageLogarithm<T>([SkipIfAllUniform, Trigger] IList<T> Uses, T result)
		where T : SettableToProduct<T>, SettableTo<T>, SettableToUniform
		{
			return ReplicateOp_NoDivide.DefAverageConditional(Uses, result);
		}
#if SpecializeArrays
		// must have upward Trigger to match the Trigger on UsesEqualDef.UsesAverageLogarithm
		/// <summary>
		/// VMP message to 'Def'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Def' with 'Uses' integrated out.
		/// The formula is <c>sum_Uses p(Uses) factor(Uses,Def,Marginal)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		[MultiplyAll]
		public static T DefAverageLogarithm<T>([SkipIfAllUniform, Trigger] T[] Uses, T result)
		where T : SettableToProduct<T>, SettableTo<T>, SettableToUniform
		{
			return ReplicateOp_NoDivide.DefAverageConditional(Uses, result);
		}
#endif
#else
			// must have upward Trigger to match the Trigger on UsesEqualDef.UsesAverageLogarithm
        public static T DefAverageLogarithm<T,TUses>([SkipIfAllUniform,Trigger] IList<TUses> Uses, T result)
            where T : SettableToProduct<TUses>, SettableTo<TUses>, TUses, SettableToUniform
        {
            return DefAverageConditional(Uses, result);
        }
#if SpecializeArrays
			// must have upward Trigger to match the Trigger on UsesEqualDef.UsesAverageLogarithm
        public static T DefAverageLogarithm<T,TUses>([SkipIfAllUniform,Trigger] TUses[] Uses, T result)
            where T : SettableToProduct<TUses>, SettableTo<TUses>, TUses, SettableToUniform
        {
            return DefAverageConditional(Uses, result);
        }
#endif
#endif
	}

	[FactorMethod(typeof(Factor), "Replicate<>")]
	[Quality(QualityBand.Mature)]
	public static class ReplicateMaxOp
	{
		/// <summary>
		/// Max product message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		public static T UsesMaxConditional<T>([AllExceptIndex] IList<T> Uses, [SkipIfUniform] T Def, int resultIndex, T result)
		where T : SettableToProduct<T>, SettableTo<T>
		{
			T res = ReplicateOp_NoDivide.UsesAverageConditional<T>(Uses, Def, resultIndex, result);
			if (res is UnnormalizedDiscrete) ((UnnormalizedDiscrete)(object)res).SetMaxToZero();
			return res;
		}
		[Skip]
		public static T UsesMaxConditionalInit<T>([IgnoreDependency] T Def, int resultIndex)
			where T : ICloneable
		{
			return (T)Def.Clone();
		}
		/// <summary>
		/// Max product message to 'Def'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Uses"/> is not a proper distribution</exception>
		public static T DefMaxConditional<T>([SkipIfAllUniform] IList<T> Uses, T result)
			where T : SettableToProduct<T>, SettableTo<T>, SettableToUniform
		{
			return ReplicateOp_NoDivide.DefAverageConditional<T>(Uses, result);
		}
		/// <summary>
		/// Max product message to 'Marginal'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		public static T MarginalMaxConditional<T>(IList<T> Uses, [SkipIfUniform] T Def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			T res = ReplicateOp_NoDivide.MarginalAverageConditional<T>(Uses, Def, result);
			if (res is UnnormalizedDiscrete) ((UnnormalizedDiscrete)(object)res).SetMaxToZero();
			return res;
		}
	}


}
