// (C) Copyright 2008 Microsoft Research Cambridge
#define FAST

using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
    /// <summary>
    /// Provides outgoing messages for <see cref="Factor.BernoulliFromBoolean(bool,double,double)"/>, given random arguments to the function.
    /// </summary>
    /// <exclude/>
    [FactorMethod(new string[] { "Sample", "Choice", "ProbTrue0", "ProbTrue1" }, typeof(Factor), "BernoulliFromBoolean", typeof(bool), typeof(double), typeof(double))]
	public static class BernoulliFromBooleanOp
	{
        /// <summary>
        /// Gibbs message to 'sample'.
        /// </summary>
        /// <param name="choice">Constant value for 'choice'.</param>
        /// <param name="probTrue0">Constant value for 'probTrue0'.</param>
        /// <param name="probTrue1">Constant value for 'probTrue1'.</param>
        /// <returns>The outgoing Gibbs message to the 'sample' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
        /// </para></remarks>
        public static Bernoulli SampleConditional(bool choice, double probTrue0, double probTrue1)
		{
            Bernoulli result = new Bernoulli();
			result.SetProbTrue(choice ? probTrue1 : probTrue0);
			return result;
		}

		/// <summary>
		/// EP message to 'sample'.
		/// </summary>
		/// <param name="choice">Constant value for 'choice'.</param>
		/// <param name="probTrue0">Constant value for 'probTrue0'.</param>
		/// <param name="probTrue1">Constant value for 'probTrue1'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'sample'.
		/// The formula is <c>int f(sample,x) q(x) dx</c> where <c>x = (choice,probTrue0,probTrue1)</c>.
		/// </para></remarks>
        public static Bernoulli SampleAverageConditional(bool choice, double probTrue0, double probTrue1)
		{
			return SampleConditional(choice, probTrue0, probTrue1);
		}

        /// <summary>
        /// VMP message to 'sample'.
        /// </summary>
        /// <param name="choice">Constant value for 'choice'.</param>
        /// <param name="probTrue0">Constant value for 'probTrue0'.</param>
        /// <param name="probTrue1">Constant value for 'probTrue1'.</param>
        /// <returns>The outgoing VMP message to the 'sample' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'sample'.
        /// The formula is <c>int log(f(sample,x)) q(x) dx</c> where <c>x = (choice,probTrue0,probTrue1)</c>.
        /// </para></remarks>
		public static Bernoulli SampleAverageLogarithm(bool choice, double probTrue0, double probTrue1)
		{
			return SampleConditional(choice, probTrue0, probTrue1);
		}

        /// <summary>
        /// Gibbs message to 'choice'.
        /// </summary>
        /// <param name="sample">Constant value for 'sample'.</param>
        /// <param name="probTrue0">Constant value for 'probTrue0'.</param>
        /// <param name="probTrue1">Constant value for 'probTrue1'.</param>
        /// <returns>The outgoing Gibbs message to the 'choice' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the factor viewed as a function of 'choice' conditioned on the given values.
        /// </para></remarks>
		public static Bernoulli ChoiceConditional(bool sample, double probTrue0, double probTrue1)
		{
            Bernoulli result = new Bernoulli();
			if (probTrue0 == 0 || probTrue1 == 0) throw new ArgumentException("probTrue is zero");
			if (sample) {
				double sum = probTrue0 + probTrue1;
				if(sum == 0.0) throw new AllZeroException();
				else result.SetProbTrue(probTrue1 / sum);
			} else {
				double sum = 2 - probTrue1 - probTrue0;
				if(sum == 0.0) throw new AllZeroException();
				else result.SetProbTrue((1 - probTrue1) / sum);
			}
			return result;
		}

        /// <summary>
        /// EP message to 'choice'.
        /// </summary>
        /// <param name="sample">Constant value for 'sample'.</param>
        /// <param name="probTrue0">Constant value for 'probTrue0'.</param>
        /// <param name="probTrue1">Constant value for 'probTrue1'.</param>
        /// <returns>The outgoing EP message to the 'choice' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'choice'.
        /// The formula is <c>int f(choice,x) q(x) dx</c> where <c>x = (sample,probTrue0,probTrue1)</c>.
        /// </para></remarks>
		public static Bernoulli ChoiceAverageConditional(bool sample, double probTrue0, double probTrue1)
		{
			return ChoiceConditional(sample,probTrue0,probTrue1);
		}

        /// <summary>
        /// VMP message to 'choice'.
        /// </summary>
        /// <param name="sample">Constant value for 'sample'.</param>
        /// <param name="probTrue0">Constant value for 'probTrue0'.</param>
        /// <param name="probTrue1">Constant value for 'probTrue1'.</param>
        /// <returns>The outgoing VMP message to the 'choice' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'choice'.
        /// The formula is <c>int log(f(choice,x)) q(x) dx</c> where <c>x = (sample,probTrue0,probTrue1)</c>.
        /// </para></remarks>
        public static Bernoulli ChoiceAverageLogarithm(bool sample, double probTrue0, double probTrue1)
		{
			return ChoiceConditional(sample,probTrue0,probTrue1);
		}

        /// <summary>
        /// EP message to 'sample'.
        /// </summary>
        /// <param name="choice">Incoming message from 'choice'.</param>
        /// <param name="probTrue0">Constant value for 'probTrue0'.</param>
        /// <param name="probTrue1">Constant value for 'probTrue1'.</param>
        /// <returns>The outgoing EP message to the 'sample' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'sample'.
        /// The formula is <c>int f(sample,x) q(x) dx</c> where <c>x = (choice,probTrue0,probTrue1)</c>.
        /// </para></remarks>
        public static Bernoulli SampleAverageConditional(Bernoulli choice, double probTrue0, double probTrue1)
		{
            Bernoulli result = new Bernoulli();
			if(choice.IsPointMass) return SampleConditional(choice.Point,probTrue0,probTrue1);
#if FAST
			result.SetProbTrue(choice.GetProbFalse() * probTrue0 + choice.GetProbTrue() * probTrue1);
#else
			// This method is more numerically stable but slower.
			// let oX = log(p(X)/(1-p(X))
			// let oY = log(p(Y)/(1-p(Y))
			// oX = log( (TT*sigma(oY) + TF*sigma(-oY))/(FT*sigma(oY) + FF*sigma(-oY)) )
			//    = log( (TT*exp(oY) + TF)/(FT*exp(oY) + FF) )
			//    = log( (exp(oY) + TF/TT)/(exp(oY) + FF/FT) ) + log(TT/FT)
			// ay = log(TF/TT)
			// by = log(FF/FT)
			// offset = log(TT/FT)
			if (probTrue0 == 0 || probTrue1 == 0) throw new ArgumentException("probTrue is zero");
			double ay = Math.Log(probTrue0 / probTrue1);
			double by = Math.Log((1 - probTrue0) / (1 - probTrue1));
			double offset = MMath.Logit(probTrue1);
			result.LogOdds = MMath.DiffLogSumExp(choice.LogOdds, ay, by) + offset;
#endif
			return result;
		}

        /// <summary>
        /// EP message to 'choice'.
        /// </summary>
        /// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
        /// <param name="probTrue0">Constant value for 'probTrue0'.</param>
        /// <param name="probTrue1">Constant value for 'probTrue1'.</param>
        /// <returns>The outgoing EP message to the 'choice' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'choice'.
        /// The formula is <c>int f(choice,x) q(x) dx</c> where <c>x = (sample,probTrue0,probTrue1)</c>.
        /// </para></remarks>
        /// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static Bernoulli ChoiceAverageConditional([SkipIfUniform] Bernoulli sample, double probTrue0, double probTrue1)
		{
            Bernoulli result = new Bernoulli();
			if(sample.IsPointMass) return ChoiceConditional(sample.Point,probTrue0,probTrue1);
#if FAST
			double p1 = sample.GetProbFalse() * (1 - probTrue1) + sample.GetProbTrue() * probTrue1;
			double p0 = sample.GetProbFalse() * (1 - probTrue0) + sample.GetProbTrue() * probTrue0;
			double sum = p0+p1;
			if(sum == 0.0) throw new AllZeroException();
			else result.SetProbTrue(p1 / sum);
#else
			// This method is more numerically stable but slower.
			// ax = log(FT/TT)
			// bx = log(FF/TF)
			// offset = log(TT/TF)
			if (probTrue0 == 0 || probTrue1 == 0) throw new ArgumentException("probTrue is zero");
			double ax = -MMath.Logit(probTrue1);
			double bx = -MMath.Logit(probTrue0);
			double offset = Math.Log(probTrue1 / probTrue0);
			result.LogOdds = MMath.DiffLogSumExp(sample.LogOdds, ax, bx) + offset;
#endif
			return result;
		}

        /// <summary>
        /// VMP message to 'sample'.
        /// </summary>
        /// <param name="choice">Incoming message from 'choice'.</param>
        /// <param name="probTrue0">Constant value for 'probTrue0'.</param>
        /// <param name="probTrue1">Constant value for 'probTrue1'.</param>
        /// <returns>The outgoing VMP message to the 'sample' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'sample'.
        /// The formula is <c>int log(f(sample,x)) q(x) dx</c> where <c>x = (choice,probTrue0,probTrue1)</c>.
        /// </para></remarks>
        public static Bernoulli SampleAverageLogarithm(Bernoulli choice, double probTrue0, double probTrue1)
		{
            Bernoulli result = new Bernoulli();
			if(choice.IsPointMass) return SampleConditional(choice.Point,probTrue0,probTrue1);
			// log(p(X=true)/p(X=false)) = sum_k p(Y=k) log(ProbTrue[k]/(1-ProbTrue[k]))
			result.LogOdds = choice.GetProbFalse() * MMath.Logit(probTrue0) + choice.GetProbTrue() * MMath.Logit(probTrue1);
			return result;
		}

        /// <summary>
        /// VMP message to 'choice'.
        /// </summary>
        /// <param name="sample">Incoming message from 'sample'.</param>
        /// <param name="probTrue0">Constant value for 'probTrue0'.</param>
        /// <param name="probTrue1">Constant value for 'probTrue1'.</param>
        /// <returns>The outgoing VMP message to the 'choice' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'choice'.
        /// The formula is <c>int log(f(choice,x)) q(x) dx</c> where <c>x = (sample,probTrue0,probTrue1)</c>.
        /// </para></remarks>
		public static Bernoulli ChoiceAverageLogarithm(Bernoulli sample, double probTrue0, double probTrue1)
		{
            Bernoulli result = new Bernoulli();
			if(sample.IsPointMass) return ChoiceConditional(sample.Point,probTrue0,probTrue1);
			// p(Y=k) =propto ProbTrue[k]^p(X=true) (1-ProbTrue[k])^p(X=false)
			// log(p(Y=true)/p(Y=false)) = p(X=true)*log(ProbTrue[1]/ProbTrue[0]) + p(X=false)*log((1-ProbTrue[1])/(1-ProbTrue[0]))
			//                           = p(X=false)*(log(ProbTrue[0]/(1-ProbTrue[0]) - log(ProbTrue[1]/(1-ProbTrue[1]))) + log(ProbTrue[1]/ProbTrue[0])
			if (probTrue0 == 0 || probTrue1 == 0) throw new ArgumentException("probTrue is zero");
			result.LogOdds = sample.GetProbTrue() * Math.Log(probTrue1 / probTrue0) + sample.GetProbFalse() * Math.Log((1 - probTrue1) / (1 - probTrue0));
			return result;
		}

	}

    /// <summary>
    /// Provides outgoing messages for <see cref="Factor.BernoulliFromBoolean(bool, double[])"/>, given random arguments to the function.
    /// </summary>
    [FactorMethod(new string[] { "Sample", "Choice", "ProbTrue" }, typeof(Factor), "BernoulliFromBoolean", typeof(bool), typeof(double[]))]
	public static class BernoulliFromBooleanArray
	{
        /// <summary>
        /// EP message to 'sample'.
        /// </summary>
        /// <param name="choice">Incoming message from 'choice'.</param>
        /// <param name="probTrue">Constant value for 'probTrue'.</param>
        /// <returns>The outgoing EP message to the 'sample' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'sample'.
        /// The formula is <c>int f(sample,x) q(x) dx</c> where <c>x = (choice,probTrue)</c>.
        /// </para></remarks>
        public static Bernoulli SampleAverageConditional(Bernoulli choice, double[] probTrue)
		{
			return BernoulliFromBooleanOp.SampleAverageConditional(choice, probTrue[0], probTrue[1]);
		}

        /// <summary>
        /// EP message to 'sample'.
        /// </summary>
        /// <param name="choice">Constant value for 'choice'.</param>
        /// <param name="probTrue">Constant value for 'probTrue'.</param>
        /// <returns>The outgoing EP message to the 'sample' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'sample'.
        /// The formula is <c>int f(sample,x) q(x) dx</c> where <c>x = (choice,probTrue)</c>.
        /// </para></remarks>
        public static Bernoulli SampleAverageConditional(bool choice, double[] probTrue)
		{
			return SampleAverageConditional(Bernoulli.PointMass(choice), probTrue);
		}

        /// <summary>
        /// EP message to 'choice'.
        /// </summary>
        /// <param name="sample">Incoming message from 'sample'.</param>
        /// <param name="probTrue">Constant value for 'probTrue'.</param>
        /// <returns>The outgoing EP message to the 'choice' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'choice'.
        /// The formula is <c>int f(choice,x) q(x) dx</c> where <c>x = (sample,probTrue)</c>.
        /// </para></remarks>
        public static Bernoulli ChoiceAverageConditional(Bernoulli sample, double[] probTrue)
		{
			return BernoulliFromBooleanOp.ChoiceAverageConditional(sample, probTrue[0], probTrue[1]);
		}

        /// <summary>
        /// VMP message to 'sample'.
        /// </summary>
        /// <param name="choice">Incoming message from 'choice'.</param>
        /// <param name="probTrue">Constant value for 'probTrue'.</param>
        /// <returns>The outgoing VMP message to the 'sample' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'sample'.
        /// The formula is <c>int log(f(sample,x)) q(x) dx</c> where <c>x = (choice,probTrue)</c>.
        /// </para></remarks>
        public static Bernoulli SampleAverageLogarithm(Bernoulli choice, double[] probTrue)
		{
			return BernoulliFromBooleanOp.SampleAverageLogarithm(choice, probTrue[0], probTrue[1]);
		}

        /// <summary>
        /// VMP message to 'choice'.
        /// </summary>
        /// <param name="sample">Incoming message from 'sample'.</param>
        /// <param name="probTrue">Constant value for 'probTrue'.</param>
        /// <returns>The outgoing VMP message to the 'choice' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'choice'.
        /// The formula is <c>int log(f(choice,x)) q(x) dx</c> where <c>x = (sample,probTrue)</c>.
        /// </para></remarks>
        public static Bernoulli ChoiceAverageLogarithm(Bernoulli sample, double[] probTrue)
		{
			return BernoulliFromBooleanOp.ChoiceAverageLogarithm(sample, probTrue[0], probTrue[1]);
		}

	}

}
