// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Power plate factor method
	/// </summary>
	[Hidden]
	public static class PowerPlate
	{
		/// <summary>
		/// Copy a value from outside to the inside of a power plate.
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="value"></param>
		/// <param name="exponent"></param>
		/// <returns>A copy of value.</returns>
		public static T Enter<T>([IsReturned] T value, double exponent) { return value; }
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="PowerPlate.Enter{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(PowerPlate), "Enter<>")]
	[Quality(QualityBand.Experimental)]
	public static class PowerPlateOp
	{
		/// <summary>
		/// EP message to 'value'.
		/// </summary>
		/// <param name="enter">Incoming message from 'enter'.</param>
		/// <param name="exponent">Constant value for 'exponent'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'value'.
		/// The formula is <c>int f(value,x) q(x) dx</c> where <c>x = (enter,exponent)</c>.
		/// </para></remarks>
		public static T ValueAverageConditional<T>([SkipIfUniform] T enter, double exponent, T result)
			where T : SettableToPower<T>
		{
			result.SetToPower(enter, exponent);
			return result;
		}

		/// <summary>
		/// EP message to 'enter'.
		/// </summary>
		/// <param name="enter">Incoming message from 'enter'.</param>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="exponent">Constant value for 'exponent'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'enter'.
		/// The formula is <c>int f(enter,x) q(x) dx</c> where <c>x = (value,exponent)</c>.
		/// </para></remarks>
		public static T EnterAverageConditional<T>(T enter, [SkipIfUniform] T value, double exponent, T result)
			where T : SettableToPower<T>, SettableToProduct<T>
		{
			if (exponent == 0) {
				// it doesn't matter what we return in this case, so we return something proper
				// to avoid spurious improper message exceptions
				result.SetToPower(value, 1.0);
			} else {
				// to_enter = value*enter^(exponent-1)
				result.SetToPower(enter, exponent - 1);
				result.SetToProduct(value, result);
			}
			return result;
		}
		[Skip]
		public static T EnterAverageConditionalInit<T>([IgnoreDependency] T value)
			where T : ICloneable
		{
			return (T)value.Clone();
		}

		public static double LogEvidenceRatio<T>(T enter, T value, double exponent, [Fresh] T to_enter)	
			where T : CanGetLogAverageOf<T>, CanGetLogAverageOfPower<T>
		{ 
			// qnot(x) =propto q(x)/m_out(x)
			// qnot2(x) =propto q(x)/m_out(x)^n
			// the interior of the plate sends (int_x qnot(x) f(x) dx)^n
			// which is missing (int_x qnot2(x) m_out(x)^n dx)/(int_x qnot(x) m_out(x) dx)^n
			// this factor sends the missing piece, where:
			// enter = m_out(x)
			// to_enter = qnot(x)
			// value = qnot2(x)
			return value.GetLogAverageOfPower(enter, exponent) - exponent*to_enter.GetLogAverageOf(enter);
		}

		/// <summary>
		/// Returns 0.
		/// </summary>
		/// <returns></returns>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'value'.
		/// </summary>
		/// <param name="enter">Incoming message from 'enter'.</param>
		/// <param name="exponent">Constant value for 'exponent'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		public static T ValueAverageLogarithm<T>([SkipIfUniform] T enter, double exponent, T result)
			where T : SettableToPower<T>
		{
			return ValueAverageConditional<T>(enter, exponent, result);
		}

		/// <summary>
		/// VMP message to 'enter'.
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		public static T EnterAverageLogarithm<T>([IsReturned] T value)
		{
			return value;
		}

	}

	/// <summary>
	/// Damp factor methods
	/// </summary>
	public static class Damp
	{
		/// <summary>
		/// Copy a value and damp the backward message.
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="value"></param>
		/// <param name="stepsize">1.0 means no damping, 0.0 is infinite damping.</param>
		/// <returns></returns>
		/// <remarks>
		/// If you use this factor, be sure to increase the number of algorithm iterations appropriately.
		/// The number of iterations should increase according to the reciprocal of stepsize.
		/// </remarks>
		public static T Backward<T>([IsReturned] T value, double stepsize) { return value; }
		/// <summary>
		/// Copy a value and damp the forward message.
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="value"></param>
		/// <param name="stepsize">1.0 means no damping, 0.0 is infinite damping.</param>
		/// <returns></returns>
		/// <remarks>
		/// If you use this factor, be sure to increase the number of algorithm iterations appropriately.
		/// The number of iterations should increase according to the reciprocal of stepsize.
		/// </remarks>
		public static T Forward<T>([IsReturned] T value, double stepsize) { return value; }
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Damp.Backward{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Damp), "Backward<>")]
	[Quality(QualityBand.Experimental)]
	public static class DampBackwardOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(backward,value,stepsize))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio()
		{
			return 0.0;
		}

		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="backward">Incoming message from 'backward'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="stepsize">Constant value for 'stepsize'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(backward) p(backward) factor(backward,value,stepsize)]/p(value)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="backward"/> is not a proper distribution</exception>
		public static Distribution ValueAverageConditional<Distribution>([SkipIfUniform] Distribution backward, double stepsize, Distribution result)
			where Distribution : SettableToPower<Distribution>, SettableToProduct<Distribution>
		{
			// damp the backward message.
			// we assume result holds the last message to value.
			// result = result^(1-stepsize) * backward^stepsize
			result.SetToPower(result, (1-stepsize)/stepsize);
			result.SetToProduct(result, backward);
			result.SetToPower(result, stepsize);
			return result;
		}

		/// <summary>
		/// EP message to 'backward'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'backward' as the random arguments are varied.
		/// The formula is <c>proj[p(backward) sum_(value) p(value) factor(backward,value,stepsize)]/p(backward)</c>.
		/// </para></remarks>
		public static Distribution BackwardAverageConditional<Distribution>([IsReturned] Distribution value)
		{
			return value;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Damp.Forward{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Damp), "Forward<>")]
	[Quality(QualityBand.Experimental)]
	public static class DampForwardOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(forward,value,stepsize))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio()
		{
			return 0.0;
		}

		/// <summary>
		/// EP message to 'forward'
		/// </summary>
		/// <param name="value">Incoming message from 'value'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="stepsize">Constant value for 'stepsize'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'forward' as the random arguments are varied.
		/// The formula is <c>proj[p(forward) sum_(value) p(value) factor(forward,value,stepsize)]/p(forward)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="value"/> is not a proper distribution</exception>
		public static Distribution ForwardAverageConditional<Distribution>([SkipIfUniform] Distribution value, double stepsize, Distribution result)
			where Distribution : SettableToPower<Distribution>, SettableToProduct<Distribution>
		{
			// damp the backward message.
			// we assume result holds the last message to value.
			// result = result^(1-stepsize) * backward^stepsize
			result.SetToPower(result, (1-stepsize)/stepsize);
			result.SetToProduct(result, value);
			result.SetToPower(result, stepsize);
			return result;
		}

		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="forward">Incoming message from 'forward'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(forward) p(forward) factor(forward,value,stepsize)]/p(value)</c>.
		/// </para></remarks>
		public static Distribution ValueAverageConditional<Distribution>([IsReturned] Distribution forward, Distribution result)
			where Distribution : SettableTo<Distribution>
		{
			result.SetTo(forward);
			return result;
		}
	}
}
