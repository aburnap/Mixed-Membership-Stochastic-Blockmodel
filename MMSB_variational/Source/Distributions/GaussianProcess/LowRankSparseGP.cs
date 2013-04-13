// LowRankSparseGP class. This is a sparseGP parameterised by
// x_i, y_i, lambda_i. This is used for message passing out of
// a GP function evaluator to a SparseGP variable
// Author: John Guiver
// (C) Copyright 2007 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Distributions.Kernels;

namespace MicrosoftResearch.Infer.Distributions.GaussianProcess
{
    /// <summary>
    /// This class is a low-rank (x,y,lambda)-parameterised sparse GP.
    /// This is used for message passing out of a GP function
    /// evaluator to a SparseGP variable
    /// </summary>
    public class LowRankSparseGP : SparseGPBase, SettableTo<LowRankSparseGP>
    {
        #region LowRankSparseGP parameters
        protected Vector xi;
        /// <summary>
        /// Input data point on which this function distribution is based
        /// </summary>
        virtual public Vector X
        {
            get
            {
                return xi;
            }

            set
            {
                xi = value;
                DoRecalculate();
            }
        }

        protected double yi;
        /// <summary>
        /// Y parameter
        /// </summary>
        virtual public double Y
        {
            get
            {
                return yi;
            }

            set
            {
                yi = value;
            }
        }

        protected double lambInv;
        /// <summary>
        /// Lambda inverse
        /// </summary>
        virtual public double LambdaInv
        {
            get
            {
                return lambInv;
            }

            set
            {
                lambInv = value;
            }
        }
        #endregion

        #region Calculated properties
        private Vector kBx;
        /// <summary>
        /// K(B,x). This is a calculated Vector maintained
        /// by the class
        /// </summary>
	    public Vector K_B_x
	    {
	      get
          {
              if (kBx == null)
                  kBx = KernelOf_X_B(X);
              return kBx;
          }
	    }
	
       private Vector pvec;
       /// <summary>
       /// p = Inv(K(B,B)) * K(B,x). This is a calculated Vector maintained
       /// by the class
        /// </summary>
	    public Vector P
	    {
	      get
          {
              if (pvec == null)
                  // Make sure we call the properties so they
                  // will recalculate if necessary
                  pvec = FixedParameters.InvKernelOf_B_B * K_B_x;
              return pvec;
          }
      }

      private double kxx;
      /// <summary>
      /// k(x)
      /// </summary>
      public double K_x_x
      {
          get
          {
              if (kxx == double.NaN)
                  // Make sure we call the properties so they
                  // will recalculate if necessary
                  kxx = Kernel.EvaluateX(X);
              return kxx;
          }
      }

      /// <summary>
      /// Function to signal recalculation of K(B,X) and p.
      /// This is called automatically if the fixed parameter
      /// class is swapped out, or if the kernel is changed, or
      /// if X changes. It should also be called by any external
      /// program modifies the kernel or other fixed parameters
      /// in place
      /// </summary>
      public override void DoRecalculate()
      {
          base.DoRecalculate();
          kBx = null;
          pvec = null;
          kxx = double.NaN;
      }
      #endregion

        #region IGaussianProcess Members
        /// <summary>
        /// Predictive Mean
        /// </summary>
        /// <param name="X">Input point</param>
        /// <returns>Predictive mean</returns>
        public double Mean(Vector X)
        {
            return 0.0;
        }

        /// <summary>
        /// Predictive Mean at a given list of points
        /// </summary>
        /// <param name="XList">List of inputs</param>
        /// <returns>Predictive mean vector</returns>
        public Vector Mean(IList<Vector> XList)
        {
            return null;
        }

        /// <summary>
        /// Predictive Variance at a given point
        /// </summary>
        /// <param name="X">Input</param>
        /// <returns>Predictive variance</returns>
        public double Variance(Vector X)
        {
            return 0.0;
        }

        /// <summary>
        /// Predictive Variance at a given list of points
        /// </summary>
        /// <param name="XList">List of inputs</param>
        /// <returns>Predictive covariance</returns>
        public PositiveDefiniteMatrix Covariance(IList<Vector> XList)
        {
            return null;
        }

        /// <summary>
        /// Predictive distribution at a given point
        /// </summary>
        /// <param name="X">Input</param>
        /// <returns>Predictive distribution</returns>
        public Gaussian Prediction(Vector X)
        {
            return null;
        }

        /// <summary>
        /// Predictive distribution at a given list of points
        /// </summary>
        /// <param name="XList">List of inputs</param>
        /// <returns>Predictive distribution</returns>
        public VectorGaussian Prediction(IList<Vector> XList)
        {
            return null;
        }
        #endregion

        #region Constructors
        /// <summary>
        /// Default constructor
        /// </summary>
        protected LowRankSparseGP()
            : base()
        {
            xi = null;
            yi = 0.0;
            lambInv = 1.0;
            kBx = null;
            pvec = null;
            kxx = double.NaN;
        }

        /// <summary>
        /// Construct a sparse GP, given basis etc
        /// </summary>
        /// <param name="spgf"></param>
        public LowRankSparseGP(SparseGPFixed spgf)
            : base(spgf)
        {
            int numBasis = spgf.NumberBasisPoints;
            xi = new Vector(numBasis);
            yi = 0.0;
            lambInv = 1.0;
            kBx = null;
            pvec = null;
            kxx = double.NaN;
        }
        #endregion Constructors

        #region ICloneable Members

        /// <summary>
        /// Clone. Note that the fixed parameters are just referenced
        /// </summary>
        /// <returns>The cloned object</returns>
        override public object Clone()
        {
            LowRankSparseGP result = new LowRankSparseGP();
            result.SetTo(this);
            return result;
        }

        #endregion

        #region SettableTo<LowRankSparseGP> Members
        /// <summary>
        /// Set one sparse GP to another. Everything is copied
        /// except the FixedParameters which are referenced.
        /// </summary>
        /// <param name="that">The sparse GP to copy</param>
        virtual public void SetTo(LowRankSparseGP that)
        {
            base.SetTo(that);
            xi.SetTo(that.xi);
            yi = that.yi;
            lambInv = that.lambInv;

            // Set the calculated values also, since these will be identical
            kBx.SetTo(that.kBx);
            pvec.SetTo(that.pvec);
        }
        #endregion
    }
}

