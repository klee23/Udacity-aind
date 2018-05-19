
from itertools import chain, combinations
from aimacode.planning import Action
from aimacode.utils import expr
from aimacode.utils import Expr

from layers import BaseActionLayer, BaseLiteralLayer, makeNoOp, make_node


class ActionLayer(BaseActionLayer):

    def _inconsistent_effects(self, actionA, actionB):
        """ Return True if an effect of one action negates an effect of the other

        See Also
        --------
        layers.ActionNode
        """
        #------ TODO: implement this function

        '''
        for x in actionA.effects:
            if x in actionB.effects:
                #print("ypu1 -> ",x," :: ", actionB.effects)
                #return True;
                print("-")
                
        for x in actionB.effects:
            if x in actionA.effects:
                #print("ypu2 -> ",x," :: ", actionA.effects)
                return True;    
        '''
        
        # print("--------",type(actionA.effects)) = <class 'frozenset'>
        # print("--------",type(actionA)) = <class 'layer.actionNode'>       
        
        for x in actionA.effects:
            # print("--------",type(x in actionA.effect)) = <class 'aimacode.utils.expr'>
            # print("--------",type(x.op)) = <class 'str'>
            # print("--------",type(x.args)) = <class 'tuple'>
            for y in actionB.effects:
                    #print("ypu3------> ",x," :: ", actionB.effects)
                    #print("x= ",x)
                    #print("x= ",expr(x))
                    #print("x.op",x.op )
                    #if(x.op!='~'):
                    #    print("x neg", Expr('~',x))
                    #print("x.arg---",x.args )
                    #print("x.arg",x.args[0] )
                    #print("y= ",y)
                    #print("y.op",y.op )
                    #print("y.arg",y.args )'''
                    if (x.op =='~' and y.op!='~'):
                        #print("ypu")
                        #print("x= ",x)
                        #print("y= ",y)
                        y_neg = Expr('~',y)
                        #print("~y= ",y_neg)
                        if y_neg.__eq__(x):
                            #print(x, " = ", y_neg)
                            return True;
                    if (y.op =='~' and x.op!='~'):
                        x_neg = Expr('~',x)
                        if x_neg.__eq__(y):
                            #print(y, " = ", x_neg)
                            return True;
                
        return False;    
        #raise NotImplementedError


    def _interference(self, actionA, actionB):
        """ Return True if the effects of either action negate the preconditions of the other 
        
        See Also
        --------
        layers.ActionNode
        """
        # TODO: implement this function
        #raise NotImplementedError
        for x in actionA.effects:
            for y in actionB.preconditions:
                if (x.op =='~' and y.op!='~'):
                    y_neg = Expr('~',y);
                    if y_neg.__eq__(x):
                        #print(x, " = ", y_neg);
                        return True;
                if (y.op =='~' and x.op!='~'):
                    x_neg = Expr('~',x);
                    if x_neg.__eq__(y):
                        #print(y, " = ", x_neg);
                        return True;
            
        for x in actionB.effects:
            for y in actionA.preconditions:
                if (x.op =='~' and y.op!='~'):
                    y_neg = Expr('~',y);
                    if y_neg.__eq__(x):
                        #print(x, " = ", y_neg);
                        return True;
                if (y.op =='~' and x.op!='~'):
                    x_neg = Expr('~',x);
                    if x_neg.__eq__(y):
                        #print(y, " = ", x_neg);
                        return True;
            
        return False;    

    def _competing_needs(self, actionA, actionB):
        """ Return True if the preconditions of the actions are all pairwise mutex in the parent layer 
        
        See Also
        --------
        layers.ActionNode
        layers.BaseLayer.parent_layer
        """
        # TODO: implement this function
        # raise NotImplementedError
        # print(self.parent_layer.is_mutex(actionA,actionB))
        # print(type(self.parent_layer.children))
        # print(self.parent_layer.children)
        for x in actionA.preconditions:
           for y in actionB.preconditions:
               # print("-----",type(x),"->",len(x.args),"->",x.args)
               #print(type(x),"->",x," :: ",type(y),"->",y)
               #print(type(x),"->",(x.__invert__)," :: ",type(y),"->",(y) )
                if (x.op =='~' and y.op!='~'):
                    y_neg = Expr('~',y);
                    if y_neg.__eq__(x):
                        #print(x, " = ", y_neg);
                        return True;
                if (y.op =='~' and x.op!='~'):
                    x_neg = Expr('~',x);
                    if x_neg.__eq__(y):
                        #print(y, " = ", x_neg);
                        return True;

            
        return False;    


class LiteralLayer(BaseLiteralLayer):

    def _inconsistent_support(self, literalA, literalB):
        """ Return True if all ways to achieve both literals are pairwise mutex in the parent layer

        See Also
        --------
        layers.BaseLayer.parent_layer
        """
        # TODO: implement this function
        # raise NotImplementedError
        
        #print("------------------------")
        #print(literalA," : literal A type ->",type(literalA)) #= <class 'aimacode.utils.expr'>
        #print(literalB," : literal B type ->",type(literalB)) #= <class 'aimacode.utils.expr'>
        
        # print(self, "self type -> ", type(self)) #= <class 'my_planning_graph.literal_layer'>
        # print(self.parent_layer, "self.parent_layer type -> ", type(self.parent_layer)) #= <class 'my_planning_graph.ActionLayer'>
        
        #print(self.parents, "type self.parents->",type(self.parents)) #= <class 'collection.defaultdict'>
        #print(self.parents[literalA], "type self.parents[literalA]->",type(self.parents[literalA]))
        #print(self.parents[literalB], "type self.parents[literalB]->",type(self.parents[literalB]))

        for x in self.parents[literalA]:
            for y in self.parents[literalB]:
               # print("x= ",x)
               # print("y= ",y)
               # print("is_mutex? = ",self.is_mutex(x,y))
               # print("parent_layer.is_mutex? = ",self.parent_layer.is_mutex(x,y))
                if not self.parent_layer.is_mutex(x,y):
                    return False;

        return True

    def _negation(self, literalA, literalB):
        """ Return True if two literals are negations of each other """
        # TODO: implement this function
        # raise NotImplementedError
        #print("----")
        #print("literal A = ",literalA," : args = ",literalA.args," : op = ", literalA.op)
        #print("literal B = ",literalB," : args = ",literalB.args," : op = ", literalB.op)

        if literalA.op == '~' and literalB.op!='~':
            temp_neg = Expr('~',literalB);
            #print("ned literalB",temp_neg);
            if temp_neg.__eq__(literalA):
                return True
        if literalA.op != '~' and literalB.op=='~':
            temp_neg = Expr('~',literalA);
            #print("ned literalA",temp_neg);
            if temp_neg.__eq__(literalB):
                return True            

        return False


class PlanningGraph:
    def __init__(self, problem, state, serialize=True, ignore_mutexes=False):
        """
        Parameters
        ----------
        problem : PlanningProblem
            An instance of the PlanningProblem class

        state : tuple(bool)
            An ordered sequence of True/False values indicating the literal value
            of the corresponding fluent in problem.state_map

        serialize : bool
            Flag indicating whether to serialize non-persistence actions. Actions
            should NOT be serialized for regression search (e.g., GraphPlan), and
            _should_ be serialized if the planning graph is being used to estimate
            a heuristic
        """
        self._serialize = serialize
        self._is_leveled = False
        self._ignore_mutexes = ignore_mutexes
        self.goal = set(problem.goal)

        # make no-op actions that persist every literal to the next layer
        no_ops = [make_node(n, no_op=True) for n in chain(*(makeNoOp(s) for s in problem.state_map))]
        self._actionNodes = no_ops + [make_node(a) for a in problem.actions_list]
        
        # initialize the planning graph by finding the literals that are in the
        # first layer and finding the actions they they should be connected to
        literals = [s if f else ~s for f, s in zip(state, problem.state_map)]
        layer = LiteralLayer(literals, ActionLayer(), self._ignore_mutexes)
        layer.update_mutexes()
        self.literal_layers = [layer]
        self.action_layers = []

    def h_levelsum(self):
        """ Calculate the level sum heuristic for the planning graph

        The level sum is the sum of the level costs of all the goal literals
        combined. The "level cost" to achieve any single goal literal is the
        level at which the literal first appears in the planning graph. Note
        that the level cost is **NOT** the minimum number of actions to
        achieve a single goal literal.
        
        For example, if Goal1 first appears in level 0 of the graph (i.e.,
        it is satisfied at the root of the planning graph) and Goal2 first
        appears in level 3, then the levelsum is 0 + 3 = 3.

        Hint: expand the graph one level at a time and accumulate the level
        cost of each goal.

        See Also
        --------
        Russell-Norvig 10.3.1 (3rd Edition)
        """
        # TODO: implement this function
        # raise NotImplementedError
        #print("-------------------")
        #print("goal = ", self.goal )
        #print("literal_layers = ", self.literal_layers, " -> ", type(self.literal_layers))
        #print("*******************")
        
        levelSum = 0;
        layer = 0;
        # create goal list independent of self.goal
        goalList = list(self.goal); 
        
        while not self._is_leveled:
            # while graph not fully explored.
            for literal in self.literal_layers[layer]:
                # check all literals in literal_layers
                for goal in goalList:
                    # check all goals against literals from above.
                    # if literal == goal...
                    if goal.__eq__(literal):
                        # increase level sum by layer amount
                        levelSum = levelSum + layer;
                        # remove goal from goalList so it is not checked twice.
                        goalList = [x for x in goalList if x is not goal];
                        # print("gial list=", goalList);
            # move to next layer...
            layer = layer + 1;
            self._extend();
            

        '''for x in range(len(self.literal_layers)):
            print(self.literal_layers[x]);
            print(type(self.literal_layers[x]))
            print(self.literal_layers[x]);
            for x in self.literal_layers[x].children:
                #print(self.literal_layers[x].children)
                print("children of literal layer = ",x);'''
                
        '''for thing in self.literal_layers:
            print("literal (thing.child)= ",thing.children);
            for kid in thing.children:
                print("literal layer children = ", kid)
                
        for idx in range(len(self.literal_layers)):
            print("self.literal_layers[",idx,"]= ", self.literal_layers[idx]);
            for thing in self.literal_layers[idx]:
                print("self.literal_layers[",idx,"] thing in = ",thing)'''

        '''print(" ")
        print("action_layers = ", self.action_layers, " -> ", type(self.action_layers))         
        for action in self.action_layers:
            print("action = ", action)
            for kid in action.children:
                print("action layer children = ",kid)
        
        print(" ")
        print("actionNodes = ", self._actionNodes);#  " -> ", type(self._actionNodes) = actionNode 
        for action in self._actionNodes:
            print("action node action = ", action, "type = ",type(action));
            print("action node action effect = ",action.effects)   
            print("action node action precon = ",action.preconditions)'''
        

        '''
        # we extend the whole graph first...
        #levelSum=0;
        while not self._is_leveled:
            self._extend();    
        for goals in self.goal:
            isGoalFound = False;
            for layer in range(len(self.literal_layers)):
                for var in self.literal_layers[layer]:
                    if goals.__eq__(var):
                        isGoalFound = True;
                        levelSum = levelSum + layer;
                        break;
                if isGoalFound:
                    break;'''
        
        return levelSum;

    def h_maxlevel(self):
        """ Calculate the max level heuristic for the planning graph

        The max level is the largest level cost of any single goal fluent.
        The "level cost" to achieve any single goal literal is the level at
        which the literal first appears in the planning graph. Note that
        the level cost is **NOT** the minimum number of actions to achieve
        a single goal literal.

        For example, if Goal1 first appears in level 1 of the graph and
        Goal2 first appears in level 3, then the levelsum is max(1, 3) = 3.

        Hint: expand the graph one level at a time until all goals are met.

        See Also
        --------
        Russell-Norvig 10.3.1 (3rd Edition)

        Notes
        -----
        WARNING: you should expect long runtimes using this heuristic with A*
        """
        # TODO: implement maxlevel heuristic
        # raise NotImplementedError
        levelSum = 0;
        layer = 0;
        # create goal list independent of self.goal
        goalList = list(self.goal); 
        
        while not self._is_leveled:
            # while graph not fully explored.
            for literal in self.literal_layers[layer]:
                # check all literals in literal_layers
                for goal in goalList:
                    # check all goals against literals from above.
                    # if literal == goal...
                    if goal.__eq__(literal):
                        # check to see if levelsum or layer is larger.
                        levelSum = max(levelSum, layer);
                        # remove goal from goalList so it is not checked twice.
                        goalList = [x for x in goalList if x is not goal];
                        # print("gial list=", goalList);
            # move to next layer...
            layer = layer + 1;
            self._extend();
            
        return levelSum;    

    def h_setlevel(self):
        """ Calculate the set level heuristic for the planning graph

        The set level of a planning graph is the first level where all goals
        appear such that no pair of goal literals are mutex in the last
        layer of the planning graph.

        Hint: expand the graph one level at a time until you find the set level

        See Also
        --------
        Russell-Norvig 10.3.1 (3rd Edition)

        Notes
        -----
        WARNING: you should expect long runtimes using this heuristic on complex problems
        """
        # TODO: implement setlevel heuristic
        # raise NotImplementedError
        setLevel = 0;
        layer = 0;
        # create goal list independent of self.goal
        goalList = list(self.goal); 
        
        while not self._is_leveled:
            if layer > 1000000:
                break;
            # while graph not fully explored.
            print("goal list subset of literl layer[",layer,"]? = ",set(goalList).issubset(self.literal_layers[layer]))
            if set(goalList).issubset(self.literal_layers[layer]):
                isGoalMutex = False;
                for goal in goalList:
                    for otherGoal in goalList:
                        print("goal = ",goal,"other goal = ",otherGoal,"mutex = ",self.literal_layers[layer].is_mutex(goal,otherGoal))
                        #print("inconsistent_support ?= ",self.literal_layers[layer]._inconsistent_support(goal,otherGoal))
                        #print("negation ?= ",self.literal_layers[layer]._negation(goal,otherGoal))
                        if self.literal_layers[layer].is_mutex(goal,otherGoal):
                            isGoalMutex = True;                                
                if not isGoalMutex:
                    setLevel = layer;
                    layer = 1000001;
                    break;

            # move to next layer...
            layer = layer + 1;
            self._extend();
            
        return setLevel; 
        
        
        

    ##############################################################################
    #                     DO NOT MODIFY CODE BELOW THIS LINE                     #
    ##############################################################################

    def fill(self, maxlevels=-1):
        """ Extend the planning graph until it is leveled, or until a specified number of
        levels have been added

        Parameters
        ----------
        maxlevels : int
            The maximum number of levels to extend before breaking the loop. (Starting with
            a negative value will never interrupt the loop.)

        Notes
        -----
        YOU SHOULD NOT THIS FUNCTION TO COMPLETE THE PROJECT, BUT IT MAY BE USEFUL FOR TESTING
        """
        while not self._is_leveled:
            if maxlevels == 0: break
            self._extend()
            maxlevels -= 1
        return self

    def _extend(self):
        """ Extend the planning graph by adding both a new action layer and a new literal layer

        The new action layer contains all actions that could be taken given the positive AND
        negative literals in the leaf nodes of the parent literal level.

        The new literal layer contains all literals that could result from taking each possible
        action in the NEW action layer. 
        """
        if self._is_leveled: return

        parent_literals = self.literal_layers[-1]
        parent_actions = parent_literals.parent_layer
        action_layer = ActionLayer(parent_actions, parent_literals, self._serialize, self._ignore_mutexes)
        literal_layer = LiteralLayer(parent_literals, action_layer, self._ignore_mutexes)

        for action in self._actionNodes:
            # actions in the parent layer are skipped because are added monotonically to planning graphs,
            # which is performed automatically in the ActionLayer and LiteralLayer constructors
            if action not in parent_actions and action.preconditions <= parent_literals:
                action_layer.add(action)
                literal_layer |= action.effects

                # add two-way edges in the graph connecting the parent layer with the new action
                parent_literals.add_outbound_edges(action, action.preconditions)
                action_layer.add_inbound_edges(action, action.preconditions)

                # # add two-way edges in the graph connecting the new literaly layer with the new action
                action_layer.add_outbound_edges(action, action.effects)
                literal_layer.add_inbound_edges(action, action.effects)

        action_layer.update_mutexes()
        literal_layer.update_mutexes()
        self.action_layers.append(action_layer)
        self.literal_layers.append(literal_layer)
        self._is_leveled = literal_layer == action_layer.parent_layer
