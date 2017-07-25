/*
Prototype/test code for simple neural net
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
  static double transferFunction(double x); // computer the mathematical transfer function between links
  static double transferFunctionDerivative(double x); // computes the derivative of the transfer funtion between links
  Neuron(unsigned numConnections, unsigned m_myIndex); // Neuron class constructor
  void setOutput(double value){m_output = value;}; // sets the output value of the neuron
  double getOutput()const {return m_output;};
  void feedForward(const Layer &prevLayer); // performs feedForward learning when passed outputs on prev layer
  void calculateGradient(double targetVals); // calulates the gradient when passed the desired values
  void calculateHiddenGradients(const Layer &nextLayer);
  void updateWeights (Layer & prevLayer); // updates the weight of the neuron based on output of the prev layer
  void setStaticParameters (double alpha, double eta) {m_alpha = alpha; m_eta = eta;}; // sets static vars

private:
  static double m_alpha; // momentum of the neuron, able to be tweaked
  static double m_eta; // learning rate of the network
  double m_output; // output value of neuron
  vector<connection> m_connections; // vector of connections in network
  double randWeight() { return (rand() /double(RAND_MAX));};
  unsigned m_myIndex; // index of the neuron in the layer
  double m_gradient; // current gradient value for neuron

};

// define Layer to be a vector of Neurons
typedef vector<Neuron> Layer;

// updates weights of the neuron based on output of the prev layer
void updateWeights (Layer & prevLayer) {


};

void Neuron::calculateHiddenGradients(const Layer &nextLayer){


};


// calculates the first oreder gradient when passed the target value
void Neuron::calculateGradient(double targetVal){

  double delta = targetVal - m_output;  // calculate derivative
  m_gradient = delta * Neuron::transferFunctionDerivative(m_output); // update m_gradient value w/ learning rate

};


// constructor for the Neuron class
Neuron::Neuron(unsigned numConnections, unsigned myIndex) {

  // append the connections vector with connection structs for each connection
  for (unsigned c = 0; c <= numConnections; ++c){
    m_connections.push_back( connection() );
    m_connections.back().weight = randWeight(); // initializes weight with a random value
    m_connections.back().weightDerivative = 0; // derivative starts at 0
  }
  m_myIndex = myIndex; // Neuron is gnostic of index
}


static double Neuron::transferFunction(double x){
  // this is the mathematical relationship that defines the feedforward operation
  // the only requirements are that its bounded at |1| and smooth
  // for this application we will use hyperbolic tangent
  return tanh(x);
};

static double Neuron::transferFunctionDerivative(double x){
  // derivative of tanh(x) is 1-tanh^2(x)
  return 1 - x * x; // appx = for the interval [1,-1]
};


// performs feedforward learning based on outputs from previous layer
void Neuron::feedForward(const Layer &prevLayer){

  // calculate the total output weight on this neuron
  double sum = 0.0;
  for (unsigned n = 0; n < prevLayer.size(); ++n){
    sum += prevLayer[n].getOutput() * prevLayer[n]::connection[myIndex].weight;
  }

  m_output = Neuron::transferFunction(sum); // pass sum to the transfer function

};




//------------------------------------ Net Class ------------------------------------------------//

class Net {
public:
  Net(vector<unsigned> &topology); // constructor that accepts a reference to a network topology
  void feedForward(const vector<double> &inputVal); // feeds the input values through the network
  void getResults(vector<double> &outputVal); // returns the results of the back propogation
  void backPropogate(const vector<double> &desiredVal); // back propogates the error through the network
  vector<Layer>::iterator layerItr(){return m_layers.begin();} // returns an iterator pointing to the start of the layer
  void setStaticParameters(double alpha, double eta); // sets the static parameters for the Neuron class

private:
  vector<Layer> m_layers; // the vector of layers that is the network
  double m_error; // total error in backpropogation
  double m_recentAverageError; // the recent average of the error
  double m_recentAverageErrorSmoothingFactor;
};

// constructor for the Net class. Accepts a vector of unsigneds defining the network topology.
Net::Net(vector<unsigned> &topology) {

  // iterate through the vector of unsigned ints to append Layers to m_layer for each layer in network
  for ( vector<unsigned>::iterator layerNum = topology.begin(); layerNum != topology.end(); ++layerNum ) {
    m_layers.push_back(Layer());
    // for each new Layer in m_layers, add as many Neurons as is specified by the topology
    for (unsigned neuronNum =0; neuronNum<= *layerNum; ++neuronNum){
      unsigned nextLayer = layerNum == topology.end()-- ? 0 : *(layerNum++);
      m_layers.back().push_back(Neuron(nextLayer, neuronNum));
      cout << "made a Neuron!" << endl; // should output every time a neuron is added to the layer
    }
  }
};


// performs the feedforward propogation.
void Net::feedForward(const vector<double> &inputVal){ // Accepts a reference to a read-only vector of doubles containing the input values.

  // check that the input vector is the correct size
  assert(inputVal.size() == m_layers[0].size() - 1);

  // loop through the first layer input neurons
  for (unsigned i = 0; i < m_layers[0].size(); ++i ){
    // set the nth neuron in the first layer to have an output equal to the nth input vector value
    m_layers[0][i]->setOutput(inputVal[i]);
  }

  // now that the inputs are set, lets feed forward the values
  for(unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
    Layer &prevLayer = m_layers[layerNum - 1];
    for(unsigned n=0; n < ( m_layers[layerNum].size()-1 ) ; ++n){
      m_layers[layerNum][n].feedForward(prevLayer);
    }
  }
};


void Net::getResults(vector<double> &outputVal){

};

void Net::setStaticParameters(double alpha, double eta){

    assert(m_layers.size() > 0);

    m_layers[0].setStaticParameters(alpha,eta);
};




void Net::backPropogate(const vector<double> &desiredVal){

  // check that the number of desired values is equal to the size of the previous layer
  assert(desiredVal.size() == m_layers.size() -1);

  //calculate overall net error from output neurons
  m_error = 0.0;
  Layer* outputLayer = m_layers.back(); // returns a reference to the last layer

  // loop through the previous layer Neurons
  for ( unsigned n= 0; n < outputLayer->size()-1; ++n ){
    double delta = *outputLayer[n].getOutput() - desiredVal[n]; // calculate the gradient
    m_error += delta * delta; // add the square to the error
  }
  m_error = sqrt(m_error / (outputLayer->size() -1)); // find the |error| of the previous layer

  m_recentAverageError = (m_recentAverageError * m_recentAverageErrorSmoothingFactor + m_error)
  / ( m_recentAverageErrorSmoothingFactor +1 ); // calulate the recent error

  // calculate output layer gradients
  for (unsigned n = 0 ; n < outputLayer->size() -1; ++n){
    outputLayer[n].calculateGradient(desiredVal[n]);
  }

  // calculate gradients on hidden layers
  for (unsigned layerNum = (m_layers.size()-1); layerNum > 0; --layerNum){
    Layer* hiddenLayer = &m_layers[layerNum]; // get a pointer to the last layer
    Layer* nextHiddenLayer = &m_layers[layerNum + 1];

    *hiddenLayer[layerNum].calculateHiddenGradients(*nextHiddenLayer);
  }

  // for all layers from outputs to first input layer, calcualte new weights
  for (unsigned layerNum = m_layers.size() -1; layerNum >0; --layerNum ) {
    Layer layer = m_layers[layerNum];
    Layer prevLayer = m_layers[layerNum +1];

    for (unsigned n = 0; n < layer.size() -1; ++n){
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

  NeuralNet.setStaticParameters(.5,.15); // set static parameters for class Neuron




  return 0;
}
