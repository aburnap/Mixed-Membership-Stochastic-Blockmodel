// This file defines the kernel factory and associated
// kernel function creator interface. Applications can
// instantiate the factory and register the Kernel Functions
// that it wants to support. The factory is needed for loading
// GP projects, as the GP itself doesn't know which set of
// kernels it is using
// Author: John Guiver
// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer;

namespace MicrosoftResearch.Infer.Distributions.Kernels
{
    /// <summary>
    /// Kernel Factory - singleton class. It maintains a list
    /// of Kernel names and their types
    /// </summary>
    public sealed class KernelFactory
    {
        private static readonly KernelFactory instance = new KernelFactory();
        private Dictionary<string, Type> creators;

        private KernelFactory()
        {
            creators = new Dictionary<string, Type>();

            // Add in the known creators
            RegisterKernelFunction(new WhiteNoise());
            RegisterKernelFunction(new SquaredExponential());
            RegisterKernelFunction(new ARD());
            RegisterKernelFunction(new LinearKernel());
            RegisterKernelFunction(new SummationKernel());
            RegisterKernelFunction(new NNKernel());
        }

        /// <summary>
        /// Registers a kernel function. The factory is primed with stock
        /// kernel functions. This function allows clients to add in custom
        /// kernel functions
        /// </summary>
        /// <param name="ikf">Type instance</param>
        public void RegisterKernelFunction(IKernelFunctionWithParams ikf)
        {
            Type kfType = ikf.GetType();
            string typeName = kfType.Name;

            if (creators.ContainsKey(typeName))
                creators[typeName] = kfType;
            else
                creators.Add(typeName, kfType);
        }

        /// <summary>
        /// Creates a kernel function by specifying the kernel function name
        /// </summary>
        /// <param name="name">Name of the kernel function</param>
        /// <returns></returns>
        public IKernelFunctionWithParams CreateKernelFunction(string name)
        {
            if (creators.ContainsKey(name))
                return Activator.CreateInstance(creators[name]) as IKernelFunctionWithParams;
            else
                return null;
        }

        /// <summary>
        /// Kernel function factory singleton instance
        /// </summary>
        public static KernelFactory Instance
        {
            get
            {
                return instance;
            }
        }
    }
}
