// Base class for SparseGP message passing
// Author: John Guiver
// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Distributions.Kernels;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;
using System.Xml.Serialization;

namespace MicrosoftResearch.Infer.Distributions
{
	/// <summary>
	/// A base class for Gaussian process distributions
	/// </summary>
	[Serializable]
	public class GaussianProcess : IGaussianProcess, Sampleable<IFunction>, IXmlSerializable
	{
		/// <summary>
		/// Mean function
		/// </summary>
		public IFunction mean;
		/// <summary>
		/// Covariance function
		/// </summary>
		public IKernelFunction kernel;

		/// <summary>
		/// Constructor
		/// </summary>
		/// <param name="mean">Mean function</param>
		/// <param name="kernel">Kernel function</param>
		[Construction("mean", "kernel")]
		public GaussianProcess(IFunction mean, IKernelFunction kernel)
		{
			this.mean = mean;
			this.kernel = kernel;
		}

		/// <summary>
		/// This base class just returns the mean function
		/// </summary>
		/// <returns></returns>
		public IFunction Sample()
		{
			return mean;
		}

		/// <summary>
		/// This base class just returns the mean function
		/// </summary>
		/// <param name="result"></param>
		/// <returns></returns>
		public IFunction Sample(IFunction result)
		{
			return Sample();
		}

		/// <summary>
		/// Mean at a given point
		/// </summary>
		/// <param name="X"></param>
		/// <returns></returns>
		public double Mean(Vector X)
		{
			return mean.Evaluate(X);
		}

		/// <summary>
		/// Mean at a given list of points
		/// </summary>
		/// <param name="XList">List of inputs</param>
		/// <returns>Predictive mean vector</returns>
		public Vector Mean(IList<Vector> XList)
		{
			int numPoints = XList.Count;
			Vector result = Vector.Zero(numPoints);
			for (int i = 0; i < numPoints; i++) {
				result[i] = Mean(XList[i]);
			}
			return result;
		}

		/// <summary>
		/// Predictive Variance at a given point
		/// </summary>
		/// <param name="X">Input</param>
		/// <returns>Predictive variance</returns>
		public double Variance(Vector X)
		{
			return kernel.EvaluateX(X);
		}

		/// <summary>
		/// Predictive covariance at a given pair of points
		/// </summary>
		/// <param name="x"></param>
		/// <param name="y"></param>
		/// <returns></returns>
		public double Covariance(Vector x, Vector y)
		{
			return kernel.EvaluateX1X2(x, y);
		}

		/// <summary>
		/// Predictive coariance at a given list of points
		/// </summary>
		/// <param name="XList">List of inputs</param>
		/// <returns>Predictive covariance</returns>
		public PositiveDefiniteMatrix Covariance(IList<Vector> XList)
		{
			int numPoints = XList.Count;
			PositiveDefiniteMatrix kXX = new PositiveDefiniteMatrix(numPoints, numPoints);
			for (int i = 0; i < numPoints; i++) {
				kXX[i, i] = kernel.EvaluateX(XList[i]);
				for (int j = i + 1; j < numPoints; j++) {
					kXX[i, j] = kernel.EvaluateX1X2(XList[i], XList[j]);
					kXX[j, i] = kXX[i, j];
				}
			}
			return kXX;
		}

		/// <summary>
		/// Predictive distribution at a given point
		/// </summary>
		/// <param name="X">Input</param>
		/// <returns>Predictive distribution</returns>
		public Gaussian Marginal(Vector X)
		{
			return new Gaussian(Mean(X), Variance(X));
		}

		/// <summary>
		/// Predictive distribution at a given list of points
		/// </summary>
		/// <param name="XList">List of inputs</param>
		/// <returns>Predictive distribution</returns>
		public VectorGaussian Joint(IList<Vector> XList)
		{
			return new VectorGaussian(Mean(XList), Covariance(XList));
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="obj"></param>
		/// <returns></returns>
		public override bool Equals(object obj)
		{
			GaussianProcess that = obj as GaussianProcess;
			if (that == null) return false;
			return (mean == that.mean) && (kernel == that.kernel);
		}

		/// <summary>
		/// 
		/// </summary>
		/// <returns></returns>
		public override int GetHashCode()
		{
			int hash = Hash.Start;
			hash = Hash.Combine(hash, mean.GetHashCode());
			hash = Hash.Combine(hash, kernel.GetHashCode());
			return hash;
		}

		/// <summary>
		/// 
		/// </summary>
		/// <returns></returns>
		public override string ToString()
		{
			return "GaussianProcess("+mean+","+kernel+")";
		}

		#region XML Serialization
		/// <summary>
		/// Parameterless constructor needed for serialization
		/// </summary>
		protected GaussianProcess() { }

		System.Xml.Schema.XmlSchema IXmlSerializable.GetSchema()
		{
			return null;
		}

		void IXmlSerializable.ReadXml(System.Xml.XmlReader reader)
		{
			reader.Read();
			reader.ReadStartElement("Mean");
			string meanTypeName = reader.Name;
			var meanSerializer = new XmlSerializer(Type.GetType(meanTypeName));
			mean = (IFunction)meanSerializer.Deserialize(reader);
			reader.ReadEndElement();
			reader.ReadStartElement("Kernel");
			string kernelTypeName = reader.Name;
			var kernelSerializer = new XmlSerializer(Type.GetType(kernelTypeName));
			kernel = (IKernelFunction)kernelSerializer.Deserialize(reader);
			reader.ReadEndElement();
			reader.ReadEndElement();
		}

		void IXmlSerializable.WriteXml(System.Xml.XmlWriter writer)
		{
			writer.WriteStartElement("Mean");
			var meanSerializer = new XmlSerializer(mean.GetType());
			meanSerializer.Serialize(writer, mean);
			writer.WriteEndElement();
			writer.WriteStartElement("Kernel");
			var kernelSerializer = new XmlSerializer(kernel.GetType());
			kernelSerializer.Serialize(writer, kernel);
			writer.WriteEndElement();
		}
		#endregion
	}
}
