// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Discrete(int,Matrix)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(new String[] { "sample", "selector", "probs" }, typeof(Factor), "Discrete", typeof(int), typeof(Matrix))]
	[Quality(QualityBand.Experimental)]
	public static class DiscreteFromDiscreteOp
	{
		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="selector">Incoming message from 'selector'.</param>
		/// <param name="probs">Constant value for 'probs'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx)</c>
		/// where <c>x = (sample,selector,probs)</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Discrete sample, Discrete selector, Matrix probs)
		{
			return Math.Log(probs.QuadraticForm(selector.GetProbs(), sample.GetProbs()));
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="selector">Incoming message from 'selector'.</param>
		/// <param name="probs">Constant value for 'probs'.</param>
		//[Skip]
		//public static double LogEvidenceRatio(Discrete sample, Discrete selector, Matrix probs) { return 0.0; }
		public static double LogEvidenceRatio(Discrete sample, Discrete selector, Matrix probs)
		{
			// use this if the rows are not normalized
			Discrete toSample = SampleAverageConditional(selector, probs, Discrete.Uniform(sample.Dimension, sample.Sparsity));
			return LogAverageFactor(sample, selector, probs)
				-toSample.GetLogAverageOf(sample);
		}

		/// <summary>
		/// EP message to 'sample'.
		/// </summary>
		/// <param name="selector">Incoming message from 'selector'.</param>
		/// <param name="probs">Constant value for 'probs'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'sample'.
		/// The formula is <c>int f(sample,x) q(x) dx</c> where <c>x = (selector,probs)</c>.
		/// </para></remarks>
		public static Discrete SampleAverageConditional(Discrete selector, Matrix probs, Discrete result)
		{
			Vector v = result.GetWorkspace();
			v.SetToProduct(selector.GetProbs(), probs);
			result.SetProbs(v);
			return result;
		}

		/// <summary>
		/// EP message to 'selector'.
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="probs">Constant value for 'probs'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'selector'.
		/// The formula is <c>int f(selector,x) q(x) dx</c> where <c>x = (sample,probs)</c>.
		/// </para></remarks>
		public static Discrete SelectorAverageConditional(Discrete sample, Matrix probs, Discrete result)
		{
			Vector v = result.GetWorkspace();
			v.SetToProduct(probs, sample.GetProbs());
			result.SetProbs(v);
			return result;
		}
	}
}
