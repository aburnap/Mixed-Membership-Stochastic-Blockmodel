// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.BernoulliFromDiscrete"/>, given random arguments to the function.
	/// </summary>
	/// <exclude/>
	[FactorMethod(new string[] { "sample", "index", "probTrue" }, typeof(Factor), "BernoulliFromDiscrete")]
	public static class BernoulliFromDiscreteOp
	{
		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="index">Incoming message from 'index'.</param>
		/// <param name="probTrue">Constant value for 'probTrue'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx)</c>
		/// where <c>x = (sample,index,probTrue)</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool sample, Discrete index, double[] probTrue)
		{
			double p = 0;
			for (int i = 0; i < index.Dimension; i++)
			{
				p += probTrue[i] * index[i];
			}
			if (!sample) p = 1-p;
			return Math.Log(p);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="index">Incoming message from 'index'.</param>
		/// <param name="probTrue">Constant value for 'probTrue'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx / int ftilde(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx / int ftilde(x) qnotf(x) dx)</c>
		/// where <c>x = (sample,index,probTrue)</c>.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool sample, Discrete index, double[] probTrue)
		{
			return LogAverageFactor(sample, index, probTrue);
		}

		/// <summary>
		/// Gibbs message to 'sample'.
		/// </summary>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <returns>The outgoing Gibbs message to the 'sample' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Bernoulli SampleConditional(int index, double[] ProbTrue)
		{
			Bernoulli result = new Bernoulli();
			result.SetProbTrue(ProbTrue[index]);
			return result;
		}

		/// <summary>
		/// EP message to 'sample'.
		/// </summary>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'sample'.
		/// The formula is <c>int f(sample,x) q(x) dx</c> where <c>x = (index,probTrue)</c>.
		/// </para></remarks>
		public static Bernoulli SampleAverageConditional(int index, double[] ProbTrue)
		{
			return SampleConditional(index, ProbTrue);
		}

		/// <summary>
		/// VMP message to 'sample'.
		/// </summary>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'sample'.
		/// The formula is <c>int log(f(sample,x)) q(x) dx</c> where <c>x = (index,probTrue)</c>.
		/// </para></remarks>
		public static Bernoulli SampleAverageLogarithm(int index, double[] ProbTrue)
		{
			return SampleConditional(index, ProbTrue);
		}

		/// <summary>
		/// Gibbs message to 'index'.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'index' conditioned on the given values.
		/// </para></remarks>
		public static Discrete IndexConditional(bool sample, double[] ProbTrue, Discrete result)
		{
			if (result == default(Discrete)) result = Discrete.Uniform(ProbTrue.Length);
			Vector prob = result.GetWorkspace();
			if (sample)
			{
				prob.SetTo(ProbTrue);
			}
			else
			{
				prob.SetTo(ProbTrue);
				prob.SetToDifference(1.0, prob);
			}
			result.SetProbs(prob);
			return result;
		}

		/// <summary>
		/// EP message to 'index'.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'index'.
		/// The formula is <c>int f(index,x) q(x) dx</c> where <c>x = (sample,probTrue)</c>.
		/// </para></remarks>
		public static Discrete IndexAverageConditional(bool sample, double[] ProbTrue, Discrete result)
		{
			return IndexConditional(sample, ProbTrue, result);
		}

		/// <summary>
		/// VMP message to 'index'.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'index'.
		/// The formula is <c>int log(f(index,x)) q(x) dx</c> where <c>x = (sample,probTrue)</c>.
		/// </para></remarks>
		public static Discrete IndexAverageLogarithm(bool sample, double[] ProbTrue, Discrete result)
		{
			return IndexConditional(sample, ProbTrue, result);
		}

		/// <summary>
		/// EP message to 'sample'.
		/// </summary>
		/// <param name="index">Incoming message from 'index'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'sample'.
		/// The formula is <c>int f(sample,x) q(x) dx</c> where <c>x = (index,probTrue)</c>.
		/// </para></remarks>
		public static Bernoulli SampleAverageConditional(Discrete index, double[] ProbTrue)
		{
			Bernoulli result = new Bernoulli();
			// E[X] = sum_Y p(Y) ProbTrue[Y]
			double p = 0;
			for (int i = 0; i < index.Dimension; i++)
			{
				p += ProbTrue[i] * index[i];
			}
			result.SetProbTrue(p);
			return result;
		}

		/// <summary>
		/// VMP message to 'sample'.
		/// </summary>
		/// <param name="index">Incoming message from 'index'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'sample'.
		/// The formula is <c>int log(f(sample,x)) q(x) dx</c> where <c>x = (index,probTrue)</c>.
		/// </para></remarks>
		public static Bernoulli SampleAverageLogarithm(Discrete index, double[] ProbTrue)
		{
			Bernoulli result = new Bernoulli();
			// E[sum_k I(Y=k) (X*log(ProbTrue[k]) + (1-X)*log(1-ProbTrue[k]))]
			// = X*(sum_k p(Y=k) log(ProbTrue[k])) + (1-X)*(sum_k p(Y=k) log(1-ProbTrue[k]))
			// p(X=true) =propto prod_k ProbTrue[k]^p(Y=k)
			// log(p(X=true)/p(X=false)) = sum_k p(Y=k) log(ProbTrue[k]/(1-ProbTrue[k]))
			double s = 0;
			for (int i = 0; i < index.Dimension; i++)
			{
				s += index[i] * MMath.Logit(ProbTrue[i]);
			}
			result.LogOdds = s;
			return result;
		}

		/// <summary>
		/// EP message to 'index'.
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'index'.
		/// The formula is <c>int f(index,x) q(x) dx</c> where <c>x = (sample,probTrue)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static Discrete IndexAverageConditional([SkipIfUniform] Bernoulli sample, double[] ProbTrue, Discrete result)
		{
			if (result == default(Discrete)) result = Discrete.Uniform(ProbTrue.Length);
			// p(Y) = ProbTrue[Y]*p(X=true) + (1-ProbTrue[Y])*p(X=false)
			Vector probs = result.GetWorkspace();
			double p = sample.GetProbTrue();
			probs.SetTo(ProbTrue);
			probs.SetToProduct(probs, 2.0 * p - 1.0);
			probs.SetToSum(probs, 1.0 - p);
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// VMP message to 'index'.
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'index'.
		/// The formula is <c>int log(f(index,x)) q(x) dx</c> where <c>x = (sample,probTrue)</c>.
		/// </para></remarks>
		public static Discrete IndexAverageLogarithm(Bernoulli sample, double[] ProbTrue, Discrete result)
		{
			if (result == default(Discrete)) result = Discrete.Uniform(ProbTrue.Length);
			// E[sum_k I(Y=k) (X*log(ProbTrue[k]) + (1-X)*log(1-ProbTrue[k]))]
			// = sum_k I(Y=k) (p(X=true)*log(ProbTrue[k]) + p(X=false)*log(1-ProbTrue[k]))
			// p(Y=k) =propto ProbTrue[k]^p(X=true) (1-ProbTrue[k])^p(X=false)
			Vector probs = result.GetWorkspace();
			double p = sample.GetProbTrue();
			probs.SetTo(ProbTrue);
			probs.SetToFunction(probs, x => Math.Pow(x, p) * Math.Pow(1.0 - x, 1.0 - p));
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// Evaluate the factor when inbox is all constants (for Gibbs sampling)
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <returns></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static double EvaluateLn(bool sample, int index, double[] ProbTrue)
		{
			double p = ProbTrue[index];
			return Math.Log(sample ? p : (1 - p));
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="index">Incoming message from 'index'.</param>
		/// <param name="ProbTrue">Constant value for 'probTrue'.</param>
		/// <returns></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static double AverageValueLn(Bernoulli sample, Discrete index, double[] ProbTrue)
		{
			double p = 0;
			for (int i = 0; i < ProbTrue.Length; i++)
			{
				p += ProbTrue[i] * index[i];
			}
			double b = sample.GetProbTrue();
			return Math.Log(b * p + (1 - b) * (1 - p));
		}
	}
}
