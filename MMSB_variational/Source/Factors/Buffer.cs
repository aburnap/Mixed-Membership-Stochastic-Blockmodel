// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;

namespace MicrosoftResearch.Infer.Factors
{
	internal static class BufferTester
	{
		[Hidden]
		public static T Copy<T>(T value) { return value; }
	}
	/// <summary>
	/// Provides outgoing messages for <see cref="BufferTester.Copy{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(BufferTester), "Copy<>")]
	[Buffers("buffer")]
	[Quality(QualityBand.Experimental)]
	public static class BufferTesterCopyOp
	{
		/// <summary>
		/// Update the buffer 'buffer'
		/// </summary>
		/// <param name="copy">Incoming message from 'copy'.</param>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static T Buffer<T>(T copy, T value, T result)
		{
			return value;
		}
		/// <summary>
		/// Initialise the buffer 'buffer'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <returns>Initial value of buffer 'buffer'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static T BufferInit<T>(T value)
		{
			return value;
		}
		/// <summary>
		/// EP message to 'copy'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="buffer">Buffer 'buffer'.</param>
		/// <returns>The outgoing EP message to the 'copy' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'copy' as the random arguments are varied.
		/// The formula is <c>proj[p(copy) sum_(value) p(value) factor(copy,value)]/p(copy)</c>.
		/// </para></remarks>
		public static T CopyAverageConditional<T>(T value, T buffer)
		{
			return value;
		}
		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="copy">Incoming message from 'copy'.</param>
		/// <returns>The outgoing EP message to the 'value' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(copy) p(copy) factor(copy,value)]/p(value)</c>.
		/// </para></remarks>
		public static T ValueAverageConditional<T>(T copy)
		{
			return copy;
		}
	}

	/// <summary>
	/// Buffer factors
	/// </summary>
	internal static class Buffer
	{
		/// <summary>
		/// Value factor
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <returns></returns>
		[Hidden]
		public static T Value<T>() { return default(T); }

		/// <summary>
		/// Infer factor
		/// </summary>
		/// <typeparam name="T"></typeparam>
		/// <param name="value"></param>
		[Hidden]
		public static void Infer<T>(T value) { }
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Buffer.Value{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Buffer), "Value<>")]
	[Quality(QualityBand.Mature)]
	internal static class BufferOp
	{
		/// <summary>
		/// EP message to 'value'.
		/// </summary>
		/// <param name="value">Incoming message from 'value'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'value' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'value'.
		/// The formula is <c>int f(value,x) q(x) dx</c> where <c>x = ()</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="value"/> is not a proper distribution</exception>
		public static T ValueAverageConditional<T>([SkipIfUniform] T value)
		{
			return value;
		}

		/// <summary>
		/// VMP message to 'value'.
		/// </summary>
		/// <param name="value">Incoming message from 'value'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'value' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'value'.
		/// The formula is <c>int log(f(value,x)) q(x) dx</c> where <c>x = ()</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="value"/> is not a proper distribution</exception>
		public static T ValueAverageLogarithm<T>([SkipIfUniform] T value)
		{
			return value;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Buffer.Infer{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Buffer), "Infer<>")]
	[Quality(QualityBand.Experimental)]
	internal static class InferOp
	{
		/// <summary>
		/// EP message to 'value'.
		/// </summary>
		/// <param name="value">Incoming message from 'value'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'value'.
		/// The formula is <c>int f(value,x) q(x) dx</c> where <c>x = (infer)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="value"/> is not a proper distribution</exception>
		public static T ValueAverageConditional<T>([SkipIfUniform] T value, T result)
	where T : SettableToUniform
		{
			result.SetToUniform();
			return result;
		}
		/// <summary>
		/// VMP message to 'value'.
		/// </summary>
		/// <param name="value">Incoming message from 'value'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'value'.
		/// The formula is <c>int log(f(value,x)) q(x) dx</c> where <c>x = (infer)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="value"/> is not a proper distribution</exception>
		public static T ValueAverageLogarithm<T>([SkipIfUniform] T value, T result)
	where T : SettableToUniform
		{
			return ValueAverageConditional(value, result);
		}
	}
}
