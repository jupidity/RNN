/*
Prototype/test code for simple neural net
based on the tutorial at https://www.youtube.com/watch?v=KkwX7FkLfug
4/4/16
*/

#include <vector>
#include <iterator>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;

// define the conneciton struct, this is so neurons will know how to pass on info
struct connection {
  double weight;
  double weightDerivative;
};

//------------------------------------ Neuron Class ---------------------------------------------//

class Neuron {
  public:
    static double transferFunction(double x);
    static double transferFunctionDerivative(double x);
    Neuron(unsigned numConnections, unsigned m_myIndex);
    void setOutput(double value){m_output = value;};
    double getOutput()const {return m_output};
    void feedForward(const Layer &prevLayer);
    void calculateGradient(double targetVals);
    void calculateHiddenGradients(const Layer &nextLayer);
private:
  static double m_alpha; // momentum of the neuron, able to be tweaked
  static double m_eta; // learning rate of the network
  double m_output;
  vector<connection> m_connections;
  double randWeight() { return (rand() /double(RAND_MAX));};
  unsigned m_myIndex;

};

Neuron::calculateHiddenGradients(const Layer &nextLayer){



}

Neuron::calculateGradient(double targetVals){

double delta = targetVal - m_output;
m_gradient = delta * Neuron::transferFunctionDerivative(m_output);

}

Neuron::Neuron(unsigned numConnections, unsigned myIndex) {
  // append the connections vector with connection structs for each connection
  for (int c = 0; c <= numConnections; ++c){
  m_connections.push_back(connection());
  m_connections.back().weight = randWeight();
  m_connections.back().weightDerivative = 0;
  }
  m_myIndex = myIndex;
}

Neuron::transferFunction(double x){
// this is the mathematical relationship that defines the feedforward operation
// the only requirements are that its bounded at |1| and smooth
return tanh(x);
}

Neuron::transferFunctionDerivative(double x){
// derivative of tanh(x) is 1-tanh^2(x)
return 1 - x * x; // appx = for the interval [1,-1]
}

Neuron::feedForward(const Layer &prevLayer){

double sum = 0.0;

for (unsigned n = 0; n < prevLayer.size(); ++n){
  sum += prevLayer[n].getOutput() * prevLayer[n]::connection[myIndex].weight;
}

m_output = Neuron::transferFunction(sum);

}

// define Layer to be a vector of Neurons
typedef vector<Neuron> Layer;


//------------------------------------ Net Class ------------------------------------------------//

class Net {
  public:
      Net(vector<unsigned> &topology); // constructor that accepts a reference to a network topology
      void feedForward(const vector<double> &inputVal); // member function that feeds the input values through the network
      void getResults(vector<double> &outputVal); // returns the results of the back propogation
      void backPropogate(const vector<double> &desiredVal); // back propogates the error through the network
      vector<Layer>::iterator layerItr(){return m_layers.begin();}
  private:
      vector<Layer> m_layers; // the vector of layers that is the network
      double m_error;
      double m_recentAverageError;
      double m_recentAverageErrorSmoothingFactor;
};

Net::Net(vector<unsigned> &topology) {
  // iterate through the vector of unsigned ints to append Layers to m_layer for each layer in network
  for ( vector<unsigned>::iterator layerNum = topology.begin(); layerNum != topology.end(); ++layerNum ) {
    m_layers.push_back(Layer());
    // for each new Layer in m_layers, add as many Neurons as is specified by the topology
    for (unsigned neuronNum =0; neuronNum<= *layerNum; ++neuronNum){
      unsigned nextLayer = layerNum == topology.end()-- ? 0 : *(layerNum++);
      m_layers.back().push_back(Neuron(nextLayer, neuronNum));
      cout << "made a Neuron!" << endl;
    }
  }
};

void Net::feedForward(const vector<double> &inputVal){
// check that the input vector is the correct size
assert(inputVal.size() == m_layers[0].size() - 1);

// loop through the first layer input neurons
for (unsigned i = 0; n < m_layers[0].size(); ++i ){
  // set the nth neuron in the first layer to have an output equal to the nth input vector value
  m_layers[0][n].setOutput(inputVal[i]);
}

// now that the inputs are set, lets feed forward the values
for(unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
  Layer &prevLayer = m_layers[layerNum - 1]);
  for(unsigned n=0; n< m_layers[layerNum].size()-1; ++n){
    m_layers.[layerNum][n].feedForward(prevLayer);
  }


}


};

void Net::getResults(vector<double> &outputVal){

};

void Net::backPropogate(const vector<double> &desiredVal){

assert(desiredVal.size() == m_layers.size() -1);

//calculate overall net error from output neurons
m_error = 0.0;
Layer &outputLayer = m_layer.back();

for (n= 0; n < outputLayer.size()-1; ++n ){
double delta = outputLayer[n].getOutput() - desiredVal[n];
m_error += delta * delta;
}
m_error = sqrt(m_error / (outputLayer.size() -1));

m_recentAverageError = (m_recentAverageError * m_recentAverageErrorSmoothingFactor + m_error)
/ ( m_recentAverageErrorSmoothingFactor +1 );

// calculate output layer gradients
for (unsigned n = 0 ; n < outputLayer.size() -1; ++n){
  outputLayer[n].calculateGradient(targetVals[n]);
}

// calculate gradients on hidden layers
for (unsigned layerNum = m_layers.size()-1; layerNum > 0; --layerNum){
  Layer &hiddenLayer = m_layers[layerNum];
  Layer &nextHiddenLayer = m_layers[layerNum + 1];

  hiddenLayer[layerNum].calculateHiddenGradients(nextHiddenLayer);
}

// for all layers from outputs to first input layer, calcualte new weights
for (unsigned LayerNum = m_layers.size() -1; layerNum >0; --layerNum ) {
  Layer &layer = m_layers[layerNum];
  Layer &prevLayer = m_layers[layerNum +1];

  for (n = 0; n < layer.size() -1; ++n){
    layer[n].updateWeights(prevLayer);
  }

}


};

int main() {

// define a network topology, here its a 4-3-1 net
vector<unsigned> topology;
topology.push_back(4);
topology.push_back(3);
topology.push_back(1);

// define the desired state variables to pass to the object
vector<double> desiredVals;
vector<double> inputVals;
vector<double> outputVals;

// create an instance of a neural net and pass the topology created above
Net NeuralNet(topology);

// adjust the alpha and eta values
//Net::Neuron::eta = .15;
//Net::Neuron::alpha = .5;

return 0;
}
