
#include "RL.h"

float EPSILON0 = 0.8;
int COLNUM = 5;
/**
constructor: initialize some v.r.
*/
RLearn::RLearn()
{
        qList.clear();
        _indexCur = 0;
        _actionSuggestIndex = 0;


}
void RLearn::rLearnInit(size_t numS, size_t numA, float ovsTarget)
{
    _numS = numS;
    _numA = numA;
    _ovsTarget = ovsTarget;
    _indexCur = 0;
    srand (time(0));
}

void RLearn::update(float state, float ovs)
{
    _action = _actionSuggest;

    // find the index of the current state and action
    cout<<" states size: "<<_states.size()<<" state: "<<_state<<endl;
    size_t stateIndexCur = findIndex(_states,_state);
    //size_t actionIndexCur = findIndex(_actions,_action);

    cout<<"* state: "<<_state<<" * action: "<<_action
    <<" * stateIndexCur: "<<stateIndexCur<<endl;

   // _indexCur = stateIndexCur*_numA + actionIndexCur;
    _indexCur = stateIndexCur;  // update the current state index

    _nextState = findDisceteValue(state, 1);
    //_action = findDisceteValue(action, 0);;
    //_qold = 0;

    // compute instant reward
    reward(ovs);

    cout<<"* _indexCur: "<<_indexCur<<" * action: "<<_action<<endl;
    cout<<"* input: "<<state<<" * nextstate: "<<_nextState<<" * Reward: "<<_reward<<endl;

    //cout<<"* Reward: "<<_reward<<endl;
}

/*void RLearn::reward(float ovs)
{
    float C1 = (ovs - _ovsTarget);
    if(fabs(C1) < 0.0005)
    {
        //_reward = 0.0001*1.0/float(pow(0.0005,2.0));
        _reward = -pow(C1*10.0,2.0)*5.0;
    }
    else if(fabs(C1) < 0.001 && C1 < 0) // small reward
        _reward = -pow(C1*10.0,2.0)*10.0;

    else if(ovs > _ovsTarget)
    {
        _reward = -pow(C1*10.0,2.0)*30.0;
    }
    else
    {
        _reward = -pow(C1*10.0,2.0)*10.0;
    }

}*/

// positive reward function
/*void RLearn::reward(float ovs)
{
    float C1 = (ovs - _ovsTarget);
    if(fabs(C1) < 0.0005)
    {
        //_reward = 0.0001*1.0/float(pow(0.0005,2.0));
        _reward = exp(-pow(C1*100.0,2.0));
    }
    else if(fabs(C1) < 0.001 && C1 < 0) // small reward
        _reward = exp(-pow(C1*100.0,2.0));

    else if(ovs > _ovsTarget)
    {
        _reward = -0.1*exp(fabs(C1*100.0));
    }
    else
    {
        _reward = exp(-pow(C1*100.0,2.0));
    }

}*/

// positive reward function
void RLearn::reward(float ovs)
{
    float C1 = (ovs - _ovsTarget);
    /*if(fabs(C1) < 0.0005)
    {
        switch(_actionSuggestIndex)
        {
            case 0:
                _reward = 2.0;
                break;
            case 1:
                _reward = 0.0;
                break;
            case 2:
                _reward = 0.0;
                break;
            case 3:
                _reward = 0.0;
                break;
            case 4:
                _reward = 0.0;
                break;
        }
    }
    else*/ if(ovs > _ovsTarget)
    {
        switch(_actionSuggestIndex)
        {
            case 0:
                _reward = 0;
                break;
            case 1:
                _reward = -10.0;
                break;
            case 2:
                _reward = -10.0;
                break;
            case 3:
                _reward = 1.0;
                break;
            case 4:
                _reward = 2.0;
                break;
        }
        //_reward = -0.1*exp(fabs(C1*100.0));
    }
    else
    {
        _reward = 0.0;
        switch(_actionSuggestIndex)
        {
            case 0:
                _reward = 0;
                break;
            case 1:
                _reward = 1.0;
                break;
            case 2:
                _reward = 2.0;
                break;
            case 3:
                _reward = -1.0;
                break;
            case 4:
                _reward = -2.0;
                break;
        }
       // _reward = exp(-pow(C1*100.0,2.0));
    }

}


void RLearn::qLearn()
{

    // Find Q(S',A') max;
    // find the index of next state
    cout<<"* q learn!"<<endl;
    size_t stateIndexNest = findIndex(_states,_nextState);

    cout<<"* stateIndexNest: "<<stateIndexNest<<endl;

    // find the action with maximum q values
    size_t index;
    QEntry qEntry = maxEntry(stateIndexNest, _nextState, index);
    _qmax = qEntry.qValue;

    cout<<"* index: "<<index<<" _actionSuggestIndex: "<<_actionSuggestIndex<<endl;


    float qcur = 0.0;
    float total = _reward + GAMMA*_qmax;

    cout<<"* _indexCur: "<<_indexCur<<" * Q value: "<<qList[_indexCur][_actionSuggestIndex].qValue<<endl;

    qcur = (1-ALPHA)*qList[_indexCur][_actionSuggestIndex].qValue+ ALPHA*total;

    // update the entry in table
    qList[_indexCur][_actionSuggestIndex].qValue = qcur;
    //_qold = qcur;

    //return qcur;
    _state = _nextState;
}

float RLearn::selectActionSuggest()
{
    cout<<"* compute EPSILON0! **********************"<<endl;
    //computeEPSILON0();
    // find the index of state
    size_t stateIndex = findIndex(_states,_state);
    //size_t stateQStart = stateIndexNext*_numA;
    size_t stateQStart = stateIndex;

    // cout<<"* stateIndex: "<<stateIndex<<" stateQStart: "<<stateQStart<<endl;
    // use epsilon greedy policy
    int randNum = (int(rand())%100);
    float randx = float(randNum)/100.0;
    cout<<"* randNum: "<<randNum<<endl;
    if(randx < EPSILON0)
    {
        // take a random action for current state
        //size_t qListSize = qList.size();
        //size_t randa = rand()%_numA;
        int randa = -1;
        while(qList[stateQStart][randa].qValue==-100 || randa == -1)
        {
            //if(stateQStart > 0 && stateQStart < _numS-1)
            randa =rand()%COLNUM;

        }

        /*else if(stateQStart == 0)   // the 1st state
            randa = rand()%COLNUM;
        else if(stateQStart == _numS-1) // the last state
            randa = 1;*/

        QEntry qEntry = qList[stateQStart][randa];
        _actionSuggest = qEntry.action;
        _actionSuggestIndex = randa;
        //_qmax = qEntry.qValue;

        //_indexCur = stateQStart+randa;

        cout<<"* randa: "<<randa<<" "<<"* randx<epsilon: "<<_actionSuggest<<endl;
    }
    else
    {
        // find the action with maximum q values
        size_t index;
        QEntry qEntry = maxEntry(stateQStart, _state, index);
        _actionSuggest = qEntry.action;
        _actionSuggestIndex = index;
        //_actionSuggestIndex =
        //_qmax = qEntry.qValue;
        //cout<<"* randx>epsilon: "<<_actionSuggest<<endl;


    }
    cout<<"* current state:++++++++++++++++++++++++++++++++++ "<<_state<<endl;
    return _actionSuggest;
}


QEntry RLearn::maxEntry(size_t beginI, float state, size_t& index)
{

    index = 0;
    QEntry qEntry;

    if(beginI >= 0)
    {

        float maxQvalue = qList[beginI][0].qValue;

        //for(size_t i = beginI; i < beginI+_numA; i++)
        for(size_t i = 0; i < COLNUM; i++)
        {
            // find all qentry with the same state
            if(qList[beginI][i].state == state && qList[beginI][i].qValue != -100)
            {
                //cout<<"* i:"<<i<<" * qList[i].qValue: "<<qList[i].qValue<<" state: "<<maxQvalue<<endl;
                if(qList[beginI][i].qValue >= maxQvalue)
                {
                    maxQvalue = qList[beginI][i].qValue;
                    index = i;
                    //cout<<"max:::"<<" "<<endl;
                }
                cout<<"* i: "<<i<<" *: "<<qList[beginI][i].qValue<<"     ";
            }

        }
    }
    else
    {
        cout<<"* RLearn::maxEntry:find nothing!"<<endl;
    }
    //_indexCur = index;

    qEntry = qList[beginI][index];
    cout<<endl<<"* max index: "<<index<<endl;

    if(index == -1)
    {
        cout<<"* RLearn::maxEntry:wrong!"<<endl;
    }

    cout<<endl;
    //cout<<"* RLearn::maxEntry:_indexCur: "<<_indexCur<<endl;

    return qEntry;


}


size_t RLearn::findIndex(vector<float>& vec, float value)
{
    size_t vecSize = vec.size();
    size_t index = -1;

    for(size_t i = 0; i < vecSize; i++)
    {
        if(vec[i] == value)
        {
            index = i;
            break;
        }
        //cout<<"vec[i]"<<vec[i]<<endl;
    }

    if(index == -1)
    {
        cout<<"* RLearn::findIndex:wrong!"<<endl;
    }

    //cout<<"* vecSize: "<<vecSize<<" * RLearn::findIndex:index+++++++++++++++++++++++++++: "<<index<<endl;
    return index;

}//end findEntry

void RLearn::initQtable(float minState, float maxState, float minAction, float maxAction
, size_t numS,  size_t numA,  float ovsTarget)
{
    //init some values first
    rLearnInit(numS,  numA,  ovsTarget);

    float granuS = (maxState- minState)/_numS;
    float granuA = (maxAction- minAction)/_numA;

    // initialize states and actions
    for(size_t i = 0; i < _numS; i++)
    {
        float stateTmp = minState + i*granuS;
        _states.push_back(stateTmp);
    }

    for(size_t i = 0; i < _numA; i++)
    {
        float actionTmp = minAction + i*granuA;
        _actions.push_back(actionTmp);
    }

    // discrete states
    vector<QEntry> qvec;

    for(size_t i = 0; i < _numS; i++)
    {
        qvec.clear();
        float stateTmp = minState + i*granuS;
        // discrete actions
        for(size_t j = 0; j < COLNUM; j++)
        {
            float actionTmp;

            // no change
            if (j == 0)
                actionTmp = stateTmp;
            else
            {
                //actionTmp = stateTmp + float(pow(-1,j)*granuA);  // indeament or decrement
                switch(j)
                {
                    case 0:
                    {
                        actionTmp = stateTmp - 0.0* granuA;
                        break;
                    }
                    case 1:
                    {
                        actionTmp = stateTmp - 2.0* granuA;
                        break;
                    }
                    case 2:
                    {
                        actionTmp = stateTmp - 1.0* granuA;
                        break;
                    }
                    case 3:
                    {
                        actionTmp = stateTmp + 1.0* granuA;
                        break;
                    }
                    case 4:
                    {
                        actionTmp = stateTmp + 2.0* granuA;
                        break;
                    }
                }
            }


            // init qList
            QEntry qEntry;
            qEntry.action = actionTmp;
            if(actionTmp < minAction || actionTmp >= maxAction)
            {
                qEntry.qValue = -100;
            }
            else
                qEntry.qValue = 0;

            if(stateTmp == 0)
            {
                if(actionTmp == 0)
                    qEntry.qValue = -100;
            }

            qEntry.state = stateTmp;

            qvec.push_back(qEntry);
            //cout<<"* state: "<<stateTmp<<" action: "<<actionTmp<<endl;
        }
        qList.push_back(qvec);
    }

    _actionSuggest = minAction;
    _nextState = _states[_states.size()/2];
    _state = _nextState;
    _indexCur = _states.size()/2;

    //cout<<" _indexCur: "<<_indexCur<<endl;
}

/**
Given a value, find the corresponding descrete one
*/
float RLearn::findDisceteValue(float value, bool isState)
{
    // find the most close one
    float minDist = 1000000;
    vector<float> vec;
    size_t vecSize = 0;
    float valueOut = value;

    if(isState == 1)    // is a state
    {
         vec = _states;
         //vecSize = _numS;

    }
    else
    {
        vec = _actions;
        //vecSize = _numA;
    }

    vecSize = vec.size();

    for(size_t i = 0; i < vecSize; i++)
    {
        float distTmp = fabs(value-vec[i]);
        if(distTmp<minDist)
        {
            valueOut = vec[i];
            minDist = distTmp;
        }
    }
    return valueOut;
}

void RLearn::computeEPSILON0()
{
    // find the index of state
    size_t stateIndex = findIndex(_states,_state);
    int nonZeroSize = 0;
    for(int i = 0; i < COLNUM; i++)
    {
        if(qList[stateIndex][i].qValue != 0 /*&& qList[stateIndex][i].qValue != -1*/)
            nonZeroSize ++;
    }
    EPSILON0 = pow(0.5,nonZeroSize);
    cout<<"* nonZeroSize: "<<nonZeroSize<<" * EPSILON0: "<<EPSILON0<<endl;
    //return EPSILON0;

}
void RLearn::clearList()
{
    for(size_t i = 0; i < _numS; i++)
    {
        for(size_t j = 0; j < COLNUM; j++)
        {
            if(qList[i][j].qValue!= -100)
            {
                qList[i][j].qValue = 0;
            }
        }
    }
}


