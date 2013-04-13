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
	/// Provides outgoing EP messages for <see cref="Factor.UsesEqualDef&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "UsesEqualDef<>")]
	[Quality(QualityBand.Mature)]
	public static class UsesEqualDefOp
	{
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
		public static double LogEvidenceRatio1<T>([SkipIfAllUniform] IList<T> Uses, T Def)
			where T : CanGetLogAverageOf<T>, SettableToProduct<T>, SettableTo<T>, ICloneable, SettableToUniform
		{
			if (Uses.Count <= 1) return 0.0;
			else {
				T toUse = (T)Def.Clone();
				T[] productBefore = new T[Uses.Count];
				T productAfter = (T)Def.Clone();
				productAfter.SetToUniform();
				double z = 0.0;
				for (int i = 0; i < Uses.Count; i++) {
					productBefore[i] = (T)Def.Clone();
					if (i > 0) productBefore[i].SetToProduct(productBefore[i-1], Uses[i-1]);
					z += productBefore[i].GetLogAverageOf(Uses[i]);
				}
				// z is now log(sum_x Def(x)*prod_i Uses[i](x))
				for (int i = Uses.Count - 1; i >= 0; i--) {
					toUse.SetToProduct(productBefore[i], productAfter);
					z -= toUse.GetLogAverageOf(Uses[i]);
					productAfter.SetToProduct(productAfter, Uses[i]);
				}
				return z;
			}
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
			if(Uses.Count <= 1) return 0.0;
			else {
				T productBefore = (T)Def.Clone();
				double z = 0.0;
				T previous_use = Def;
				for (int i = 0; i < Uses.Count; i++) {
					if (i > 0) productBefore.SetToProduct(productBefore, previous_use);
					T use = Uses[i];
					z += productBefore.GetLogAverageOf(use);
					// z is now log(sum_x Def(x)*prod_i Uses[i](x))
					z -= to_Uses[i].GetLogAverageOf(use);
					previous_use = use;
				}
				return z;
			}
		}

		/// <summary>
		/// EP message to 'Marginal'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Marginal' as the random arguments are varied.
		/// The formula is <c>proj[p(Marginal) sum_(Uses,Def) p(Uses,Def) factor(Uses,Def,Marginal)]/p(Marginal)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		// TM: SkipIfUniform on Def added as a stronger constraint, to prevent improper messages in EP.
		//[SkipIfAllUniform]
		[MultiplyAll]
		public static T MarginalAverageConditional<T>(IList<T> Uses, [SkipIfAllUniform] T Def, T result)
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
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Marginal' as the random arguments are varied.
		/// The formula is <c>proj[p(Marginal) sum_(Uses,Def) p(Uses,Def) factor(Uses,Def,Marginal)]/p(Marginal)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		// TM: SkipIfUniform on Def added as a stronger constraint, to prevent improper messages in EP.
		//[SkipIfAllUniform]
		[MultiplyAll]
		public static T MarginalAverageConditional<T>(T[] Uses, [SkipIfAllUniform] T Def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			result.SetTo(Def);
			return Distribution.SetToProductWithAll(result, Uses);
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
		// TM: SkipIfUniform on Def added as a stronger constraint, to prevent improper messages in EP.
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
		// TM: SkipIfUniform on Def added as a stronger constraint, to prevent improper messages in EP.
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
#else
			public static T DefAverageConditional<T,TUses>([SkipIfAllUniform] IList<TUses> Uses, T result)
            where T : SettableToProduct<TUses>, SettableTo<TUses>, TUses, SettableToUniform
        {
            return Distribution.SetToProductOfAll(result, Uses);
        }
#endif
	}

	/// <summary>
	/// Provides outgoing Gibbs messages for <see cref="Factor.UsesEqualDef&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "UsesEqualDef<>")]
	[Quality(QualityBand.Mature)]
	public static class UsesEqualDefGibbsOp<T>
	{
		/// <summary>
		/// Evidence message for Gibbs
		/// </summary>
		public static double GibbsEvidence<TDist>(IList<TDist> Uses, TDist Def, GibbsMarginal<TDist,T> to_marginal)
			where TDist : IDistribution<T>, Sampleable<T>, CanGetLogAverageOf<TDist>, SettableTo<TDist>, SettableToProduct<TDist>
		{
			if (Uses.Count == 1) {
				// the total evidence contribution of this variable should be Def.GetLogAverageOf(Uses[0]).
				// but since this variable is sending a sample to Def and Use, and those factors will send their own evidence contribution,
				// we need to cancel the contribution of those factors here.
				return Def.GetLogAverageOf(Uses[0]) -Def.GetLogProb(to_marginal.LastSample) -Uses[0].GetLogProb(to_marginal.LastSample);
			} else {
				//throw new ApplicationException("Gibbs Sampling does not support variables defined within a gate");
				double z = 0.0;
				TDist productBefore = (TDist)Def.Clone();
				TDist product = (TDist)Def.Clone();
				for (int i = 0; i < Uses.Count; i++) {
					if (i > 0) product.SetToProduct(productBefore, Uses[i-1]);
					z += product.GetLogAverageOf(Uses[i]);
					productBefore.SetTo(product);
				}
				// z is now log(sum_x Def(x)*prod_i Uses[i](x)), which is the desired total evidence.
				// but we must also cancel the contribution of the parent and child factors that received a sample from us.
				z -= Def.GetLogProb(to_marginal.LastSample);
				for (int i = 0; i < Uses.Count; i++) {
					z -= Uses[i].GetLogProb(to_marginal.LastSample);
				}
				return z;
			}
		}

		/// <summary>
		/// Gibbs message to 'Marginal'.
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_marginal">Previous outgoing message to 'marginal'.</param>
		/// <returns></returns>
		/// <remarks><para>
		/// The outgoing message is the product of 'Def' and 'Uses' messages.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		[Stochastic]
		//[SkipIfAllUniform("Uses","Def")]
		public static GibbsMarginal<TDist,T> MarginalGibbs<TDist>(
			IList<TDist> Uses,
			[SkipIfAnyUniform]TDist Def,
			GibbsMarginal<TDist,T> to_marginal) // must not be called 'result', because its value is used
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToRatio<TDist>, SettableTo<TDist>, Sampleable<T>
		{
			GibbsMarginal<TDist,T> result = to_marginal;
			TDist marginal = result.LastConditional;
			marginal.SetTo(Def);
			marginal = Distribution.SetToProductWithAll(marginal, Uses);
			result.LastConditional = marginal;
			// Allow a sample to be drawn from the last conditional, and add it to the sample
			// list and conditional list
			result.PostUpdate();
			return result;
		}

		[Stochastic]
		public static GibbsMarginal<TDist, T> MarginalGibbs<TDist>(
			IList<TDist> Uses,
			T Def,
			GibbsMarginal<TDist, T> to_marginal) // must not be called 'result', because its value is used
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToRatio<TDist>, SettableTo<TDist>, Sampleable<T>
		{
			GibbsMarginal<TDist,T> result = to_marginal;
			TDist marginal = result.LastConditional;
			marginal.Point = Def;
			result.LastConditional = marginal;
			// Allow a sample to be drawn from the last conditional, and add it to the sample
			// list and conditional list
			result.PostUpdate();
			return result;
		}

		[Stochastic]
		public static GibbsMarginal<TDist, T> MarginalGibbs<TDist>(
			IList<T> Uses,
			[SkipIfUniform]TDist Def,
			GibbsMarginal<TDist, T> to_marginal) // must not be called 'result', because its value is used
			where TDist : IDistribution<T>, SettableToProduct<TDist>, SettableToRatio<TDist>, SettableTo<TDist>, Sampleable<T>
		{
			if (Uses.Count > 1) throw new ArgumentException("Uses.Count > 1");
			GibbsMarginal<TDist,T> result = to_marginal;
			TDist marginal = result.LastConditional;
			marginal.Point = Uses[0];
			result.LastConditional = marginal;
			// Allow a sample to be drawn from the last conditional, and add it to the sample
			// list and conditional list
			result.PostUpdate();
			return result;
		}

		[Skip]
		public static GibbsMarginal<TDist, T> MarginalGibbsInit<TDist>([IgnoreDependency] TDist Def)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			return new GibbsMarginal<TDist, T>(Def, 100, 1, true, true, true);
		}

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
		public static T UsesGibbs<TDist>([SkipIfUniform] GibbsMarginal<TDist, T> to_marginal, int resultIndex, T result)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			return to_marginal.LastSample;
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
		public static ArrayType UsesGibbsInit<ArrayType,T>([IgnoreDependency] T Def, int count, IArrayFactory<T, ArrayType> factory)
			 where T : ICloneable
		{
			return factory.CreateArray(count, i => (T)Def.Clone());
		}

		/// <summary>
		/// Gibbs sample message to 'Def'
		/// </summary>
		/// <param name="to_marginal">Incoming message from 'Marginal'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the current Gibbs sample.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="to_marginal"/> is not a proper distribution</exception>
		public static T DefGibbs<TDist>([SkipIfUniform] GibbsMarginal<TDist, T> to_marginal, T result)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			return to_marginal.LastSample;
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
	/// Provides outgoing Gibbs messages for <see cref="Factor.UsesEqualDef&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "UsesEqualDefGibbs<>")]
	[Buffers("sample","conditional","marginalEstimator","sampleAcc","conditionalAcc")]
	[Quality(QualityBand.Mature)]
	public static class UsesEqualDefGibbsOp2<T>
	{
		[Skip]
		public static TDist ConditionalInit<TDist>([IgnoreDependency] TDist def)
			where TDist : ICloneable
		{
			return (TDist)def.Clone();
		}
		public static TDist Conditional<TDist>(IList<TDist> Uses, [SkipIfAnyUniform]TDist Def, TDist result)
			where TDist : SettableTo<TDist>, SettableToProduct<TDist>
		{
			result.SetTo(Def);
			result = Distribution.SetToProductWithAll(result, Uses);
			return result;
		}
		[Stochastic]
		public static T Sample<TDist>([IgnoreDependency] TDist def, [Proper]TDist conditional)
			where TDist : IDistribution<T>, Sampleable<T>
		{
			return conditional.Sample();
		}

		public static BurnInAccumulator<TDist> MarginalEstimatorInit<TDist,T>([IgnoreDependency] TDist to_marginal, int burnIn)
			where TDist : IDistribution<T>
		{
			Accumulator<TDist> est = (Accumulator<TDist>)ArrayEstimator.CreateEstimator<TDist,T>(to_marginal, true);
			return new BurnInAccumulator<TDist>(burnIn, 1, est);
		}
		public static TAcc MarginalEstimator<TDist,TAcc>([Proper]TDist conditional, TAcc marginalEstimator)
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
		/// Evidence message for Gibbs
		/// </summary>
		public static double GibbsEvidence<TDist>(IList<TDist> Uses, TDist Def, T sample)
			where TDist : IDistribution<T>, Sampleable<T>, CanGetLogAverageOf<TDist>, SettableTo<TDist>, SettableToProduct<TDist>
		{
			if (Uses.Count == 1) {
				// the total evidence contribution of this variable should be Def.GetLogAverageOf(Uses[0]).
				// but since this variable is sending a sample to Def and Use, and those factors will send their own evidence contribution,
				// we need to cancel the contribution of those factors here.
				return Def.GetLogAverageOf(Uses[0]) -Def.GetLogProb(sample) -Uses[0].GetLogProb(sample);
			} else {
				//throw new ApplicationException("Gibbs Sampling does not support variables defined within a gate");
				double z = 0.0;
				TDist productBefore = (TDist)Def.Clone();
				TDist product = (TDist)Def.Clone();
				for (int i = 0; i < Uses.Count; i++) {
					if (i > 0) product.SetToProduct(productBefore, Uses[i-1]);
					z += product.GetLogAverageOf(Uses[i]);
					productBefore.SetTo(product);
				}
				// z is now log(sum_x Def(x)*prod_i Uses[i](x)), which is the desired total evidence.
				// but we must also cancel the contribution of the parent and child factors that received a sample from us.
				z -= Def.GetLogProb(sample);
				for (int i = 0; i < Uses.Count; i++) {
					z -= Uses[i].GetLogProb(sample);
				}
				return z;
			}
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
		/// <exception cref="ImproperMessageException"><paramref name="marginal"/> is not a proper distribution</exception>
		public static T UsesGibbs<TDist>(TDist def, T sample, int resultIndex, T result)
			where TDist : IDistribution<T>
		{
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

		/// <summary>
		/// Gibbs sample message to 'Def'
		/// </summary>
		/// <param name="sample">Current sample.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the current Gibbs sample.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="marginal"/> is not a proper distribution</exception>
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
	/// Provides outgoing max product messages for <see cref="Factor.UsesEqualDef&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "UsesEqualDef<>")]
	[Quality(QualityBand.Mature)]
	public static class UsesEqualDefMaxOp
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
			T res = UsesEqualDefOp.UsesAverageConditional<T>(Uses, Def, resultIndex, result);
			if (res is UnnormalizedDiscrete) ((UnnormalizedDiscrete)(object)res).SetMaxToZero();
			return res;
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
			return UsesEqualDefOp.DefAverageConditional<T>(Uses, result);
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
			T res = UsesEqualDefOp.MarginalAverageConditional<T>(Uses, Def, result);
			if (res is UnnormalizedDiscrete) ((UnnormalizedDiscrete)(object)res).SetMaxToZero();
			return res;
		}
	}

	/// <summary>
	/// Provides outgoing VMP messages for <see cref="Factor.UsesEqualDef&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "UsesEqualDef<>", Default = true)]
	[Quality(QualityBand.Mature)]
	public static class UsesEqualDefVmpBufferOp
	{
		/// <summary>
		/// VMP message to 'Marginal'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'Marginal'.
		/// The formula is <c>exp(sum_(Uses,Def) p(Uses,Def) log(factor(Uses,Def,Marginal)))</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		[MultiplyAll]
		public static T MarginalAverageLogarithm<T>(IList<T> Uses, T Def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			result.SetTo(Def);
			return Distribution.SetToProductWithAll(result, Uses);
		}

		/// <summary>
		/// VMP message to 'Uses'
		/// </summary>
		/// <param name="to_marginal">Previous outgoing message to 'marginal'.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Uses' conditioned on the given values.
		/// </para></remarks>
		public static T UsesAverageLogarithm<T>([IsReturned] T to_marginal, int resultIndex, T result)
			where T : SettableTo<T>
		{
			result.SetTo(to_marginal);
			return result;
		}

		/// <summary>
		/// VMP message to 'Def'
		/// </summary>
		/// <param name="to_marginal">Previous outgoing message to 'marginal'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Def' conditioned on the given values.
		/// </para></remarks>
		public static T DefAverageLogarithm<T>([IsReturned] T to_marginal, T result)
			where T : SettableTo<T>
		{
			result.SetTo(to_marginal);
			return result;
		}
	}

	/// <summary>
	/// Provides outgoing VMP messages for <see cref="Factor.UsesEqualDef&lt;T&gt;"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "UsesEqualDef<>")]
	[Quality(QualityBand.Mature)]
	public static class UsesEqualDefVmpOp
	{
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="to_marginal">Outgoing message to 'marginal'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Uses,Def,Marginal))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor<T>([Fresh] T to_marginal)
			where T : CanGetAverageLog<T>
		{
			return -to_marginal.GetAverageLog(to_marginal);
		}

#if MinimalGenericTypeParameters
		/// <summary>
		/// VMP message to 'Marginal'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'Marginal'.
		/// The formula is <c>exp(sum_(Uses,Def) p(Uses,Def) log(factor(Uses,Def,Marginal)))</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		public static T MarginalAverageLogarithm<T>(IList<T> Uses, T Def, T result)
				where T : SettableToProduct<T>, SettableTo<T>
		{
			return UsesAverageLogarithm(Uses, Def, 0, result);
		}

		/// <summary>
		/// VMP message to 'Uses'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="resultIndex">Index of the 'Uses' array for which a message is desired</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'Uses'.
		/// The formula is <c>exp(sum_(Def) p(Def) log(factor(Uses,Def,Marginal)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		// TM: SkipIfUniform on Def added as a stronger constraint, to prevent improper messages.
		//[SkipIfAllUniform]
		[MultiplyAll]
		public static T UsesAverageLogarithm<T>(IList<T> Uses, [Proper] T Def, int resultIndex, T result)
					where T : SettableToProduct<T>, SettableTo<T>
		{
			result.SetTo(Def);
			return Distribution.SetToProductWithAll(result, Uses);
		}

		/// <summary>
		/// VMP message to 'Def'
		/// </summary>
		/// <param name="Uses">Incoming message from 'Uses'.</param>
		/// <param name="Def">Incoming message from 'Def'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'Def'.
		/// The formula is <c>exp(sum_(Uses) p(Uses) log(factor(Uses,Def,Marginal)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Def"/> is not a proper distribution</exception>
		// TM: Proper added on Def to avoid improper messages.
		//[SkipIfAllUniform]
		public static T DefAverageLogarithm<T>(IList<T> Uses, [Proper] T Def, T result)
					where T : SettableToProduct<T>, SettableTo<T>
		{
			return UsesAverageLogarithm(Uses, Def, 0, result);
		}
#else
		[SkipIfAllUniform]
		public static T MarginalAverageLogarithm<T, TUses, TDef>(IList<TUses> Uses, TDef Def, T result)
				where T : SettableToProduct<TUses>, SettableTo<TDef>, TUses
		{
			return UsesAverageLogarithm<T, TUses, TDef>(Uses, Def, 0, result);
		}

		[SkipIfAllUniform]
		public static T UsesAverageLogarithm<T, TUses, TDef>([MatchingIndexTrigger] IList<TUses> Uses, TDef Def, int resultIndex, T result)
					where T : SettableToProduct<TUses>, SettableTo<TDef>, TUses
		{
			result.SetTo(Def);
			return Distribution.SetToProductWithAll(result, Uses);
		}

		[SkipIfAllUniform]
		public static T DefAverageLogarithm<T, TUses, TDef>(IList<TUses> Uses, [Trigger] TDef Def, T result)
					where T : SettableToProduct<TUses>, SettableTo<TDef>, TUses
		{
			return UsesAverageLogarithm<T, TUses, TDef>(Uses, Def, 0, result);
		}
#endif
	}
}
