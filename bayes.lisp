;;;;;; The TABLE-BASED BAYES FILTER and THE PARTICLE FILTER
;;;;;; These are two simple filter examples.  The example world model
;;;;;; we have provided is a discrete model, which is a little silly
;;;;;; for the particle filter, and so the total amount of code for
;;;;;; the praticle filter is larger than for the discrete table-based
;;;;;; Bayes filter.  But in a continuous example, where we'd have to
;;;;;; discretize the model, the Bayes filter would get much uglier and
;;;;;; the Particle filter would stay very pretty.



;;;;;;; What you need to do:
;;;;;;;
;;;;;;; 1. Implement the Bayes and Particle Filters.  You may use different
;;;;;;;    internal functions than the templates that I have provided, but please
;;;;;;;    use the same BAYES-FILTER and PARTICLE-FILTER function forms.
;;;;;;;
;;;;;;; 2. Demonstrate that the filters produce the same results.  What happens
;;;;;;;    when you change the number of particles the particle filter uses?
;;;;;;;
;;;;;;; 3. Implement a fun new problem domain for one or both of the filters.  If you wish to
;;;;;;;    do a problem domain for FORWARD function rather than one for a filter,
;;;;;;;    you merely need to specify a single action!
;;;;;;;
;;;;;;; 4. Explore different situations in your new problem domain.


;;;;;;; OUR EXAMPLE MODEL
;;;;;;;
;;;;;;; This simple model is discrete, to keep things simple -- but in fact
;;;;;;; our particle filter should be able to be modified for a continuous
;;;;;;; world fairly trivially. 
;;;;;;;
;;;;;;; In this model we have five rooms:  0 1 2 3 4
;;;;;;; the rooms are organized in a torus (a loop).
;;;;;;;
;;;;;;; The robot can perform two actions:   :forward    and    :backward
;;;;;;; The robot can sense two possible sensations in a room:   :even   and   :odd
;;;;;;;
;;;;;;; If a robot is in room x, and tries to go forward, with 1/2 probability
;;;;;;; he will wind up in room x+1 and with 1/2 probability he will stay in room x
;;;;;;;
;;;;;;; If a robot is in an even room (0, 2, or 4) with 95% probability he will
;;;;;;; sense ":even" and with 5% probability he will sense ":odd".  The opposite
;;;;;;; occurs if he's in an odd room (1 or 3). 
;;;;;;;
;;;;;;; Hint: instead of directly accessing *action-model*, *sensor-model*, and *states*
;;;;;;; in your code later on, you might use the states, action-probability, and
;;;;;;; sensor-probability functions instead.  That way if you change your model
;;;;;;; you just need to resubmit those three functions to




;;; Parts 1 & 2:
;;;
;;; we implemented the bayes and particle filters using the provided function. On the base example runs, we got the following results:
;;;
;;; Example 1:
;;;   Bayes Filter Results: (0.7326335 0.040256143 0.19678028 0.0089635225 0.021366648)
;;;   Particle Filter Results: (0.7284 0.04        0.1999     0.0086       0.0231)
;;;
;;; Example 2:
;;;   Bayes Filter Results: (0.045372415 0.004325214 2.1597654e-4 0.8596945 0.09039186)
;;;   Particle Filter Results: (0.0427   0.0043      0.0          0.8614    0.0916)
;;;
;;; Where N, the number of particles, is equal to 10^4. One sees that with rounding, the filters (roughly) agree up to the second decimal place.
;;;
;;; Increasing N = 5*10^4 yields:
;;;   Example 1: (0.7285 0.04178 0.1992 0.00944 0.02108)
;;;   Example 2: (0.04452 0.00516 2.4e-4 0.8601 0.08998)
;;;
;;; N = 10^5:
;;;   Example 1: (0.72704726 0.040590405 0.20100202 0.00906009 0.022300223)
;;;   Example 2: (0.04439 0.00427 2.9e-4 0.86001 0.09104)
;;;
;;;
;;; which are in pretty good agreement. 
;;;
;;; We can see that for state=2 in example 2, probability starts accruing with higher N where before there was nothing. With high enough N, I expect the particle filter and the bayesian filter results to overlap completely, up to floating point error.
;;;
;;; To get these results, i.e. our answers for 1 and 2, uncomment the example model state and run the example-1-particle, example-2-bayes, etc., functions as usual.
;;;
;;; Part 3
;;;
;;; For a new problem domain, we chose to make a continuous problem setting for the particle filter.
;;; The goal was to extend the loop model from the discretized 5 rooms to allowing the robot to be
;;; anywere on the loop.
;;;
;;; To add some non-linearity, one point in the loop is blocked off by a wall. If the robot tries
;;; to pass through the wall, it bounces off elastically. 
;;; Thus the state can be given by
;;;
;;;    theta in the interval [0, 2 pi)
;;;
;;; where the theat=0 point is the wall. We allow the robot to think it can be directly in the wall
;;; to make the math easier.
;;;
;;; The sensor model:
;;;
;;; is given by "distance to the wall" - the robot knows, (with some noise) the shortest arc
;;; between it and the wall - perhaps it has invented a wonderful curved laser. We choose the sensor
;;; explicitly for the symmetry it creates - there are two points equidistant from the wall on the
;;; loop. We expect experiments to keep this symmetry.
;;;
;;; Explicitly, a sensor measurement yields
;;;
;;;    phi = (theta to wall) + epsilon
;;;
;;; where epsilon is sampled from a normal distribution with known (or guessed) standard deviation.
;;; 
;;; From phi, we can generate a weight function on theta to be a gaussian centered on phi.
;;;
;;; The action model:
;;;
;;; The robot can go forwards or backwards. Since the loop is symmetric, this has to mean something
;;; to the robot in particular. We define fowards as increasing theta and backwards as decreasing
;;; theta.
;;;
;;; If forwards, move by some average angle forwards, with some normal error (similar to sensor error).
;;; if backwards, move by perhaps a different angle backwards (maybe the motor strength is different
;;; for reversing) and add error similarly.
;;;
;;; Collisions are handled as follows. If forward would push the robot to a mean position of 2pi + k,
;;; when sampling the exact state, map any result on (2pi + k) to (2pi - k).  This means the
;;; resulting probability distribution will *not* be symmetric around the mean.
;;;
;;; Similarly for backwards, map any sample on (-k, 0) to (0,k). This will also introduce an asymmetry.
;;;
;;; The functions defining this model are outlined in the NEW MODEL section - the discrete loop is left
;;; as is in the OLD MODEL section so results from Part 1 and 2 can be reproduced. The methods to toggle
;;; between the discrete and continue model are found just above the particle filter method.
;;;
;;; Part 4:
;;;
;;; We fix N=10000.
;;;
;;; The first test we implemented was for an initial belief concentrated on the interval [pi-0.2, pi+0.2],
;;; simulating a robot driving taking the forward and backward action 5 times in succession.
;;;
;;; This interval is just about as far as we can get from the wall, so this test is mostly a sanity
;;; check that the model is correct and works well with the filter. 
;;;
;;; One (forward, backwards) pair moves the robot ~0.157 radians, since the average forward movement
;;; slightly exceeds the average backwards movement. We placed the sensor readings at the average
;;; movement from a forward/backward action respectively.
;;;
;;; Run this test with the command (example-1-continuous), with the proper model active.
;;;
;;; A sample output:
;;;
;;; Probability mass on [3.144,4.589]
;;;
;;; 3.144 to 3.216: 0.00030003
;;; 3.216 to 3.288: 0.00010001
;;; 3.288 to 3.360: 0.00050005
;;; 3.360 to 3.433: 0.00220022
;;; 3.433 to 3.505: 0.00690069
;;; 3.505 to 3.577: 0.01790179
;;; 3.577 to 3.649: 0.04070407
;;; 3.649 to 3.722: 0.06810681
;;; 3.722 to 3.794: 0.104810484
;;; 3.794 to 3.866: 0.13381338
;;; 3.866 to 3.939: 0.1530153
;;; 3.939 to 4.011: 0.1441144
;;; 4.011 to 4.083: 0.12451245
;;; 4.083 to 4.155: 0.09620962
;;; 4.155 to 4.228: 0.05480548
;;; 4.228 to 4.300: 0.03160316
;;; 4.300 to 4.372: 0.01340134
;;; 4.372 to 4.445: 0.00530053
;;; 4.445 to 4.517: 0.00120012
;;; 4.517 to 4.589: 0.00050005
;;;
;;; We expect the robot to be in position ~ (3.14+5*0.157) at the end, or around 3.927. The distribution
;;; consistently peaks ~ 3.9, which means our model passes the first test.
;;;
;;; The second test is implemented to show how powerful sensor signals around (pi) are. Even if noise
;;; in the sensor vanished, a given measurement would still leave two candidate states on the loop,
;;; with the exception of a sensor reading of pi. In a way, you get the *most* information for
;;; measurements around pi.
;;;
;;; We start by initializing an initial belief concentrated on [2 pi/3 , 3 pi/4] and [pi/2 pi/3], but
;;; run the *same* action/sensor sequence as in test 1. i.e. It is likely the robot's initial belief
;;; is quite flawed. We see if 5 (forward backward) sensor signals are enough to correct this flawed initial belief.
;;;
;;; Sample output:
;;; 
;;; 3.133 to 3.206: 0.00010001
;;; 3.206 to 3.279: 0.00020002
;;; 3.279 to 3.352: 0.00110011
;;; 3.352 to 3.425: 0.0033003301
;;; 3.425 to 3.498: 0.010001
;;; 3.498 to 3.571: 0.02370237
;;; 3.571 to 3.644: 0.046304632
;;; 3.644 to 3.717: 0.07560756
;;; 3.717 to 3.790: 0.1140114
;;; 3.790 to 3.864: 0.14571457
;;; 3.864 to 3.937: 0.15541553
;;; 3.937 to 4.010: 0.1491149
;;; 4.010 to 4.083: 0.10771077
;;; 4.083 to 4.156: 0.08320832
;;; 4.156 to 4.229: 0.04460446
;;; 4.229 to 4.302: 0.02460246
;;; 4.302 to 4.375: 0.00990099
;;; 4.375 to 4.448: 0.00370037
;;; 4.448 to 4.522: 0.0015001501
;;; 4.522 to 4.595: 0.00020002
;;;
;;; with 5 (foward backward) pairs, this distribution is nearly identical to the test-1 distribution.
;;; What if we restrict this to just 1 (forward backward) pair (commenting out the larger action-sequence)?
;;;
;;; Sample output:
;;;
;;; Probability mass on [1.602,3.221]

;;; 1.602 to 1.683: 0.000100120145
;;; 1.683 to 1.764: 0.0
;;; 1.764 to 1.845: 0.0
;;; 1.845 to 1.926: 0.0006007209
;;; 1.926 to 2.007: 0.002803364
;;; 2.007 to 2.088: 0.0073087704
;;; 2.088 to 2.169: 0.019122947
;;; 2.169 to 2.250: 0.045554664
;;; 2.250 to 2.331: 0.075390466
;;; 2.331 to 2.412: 0.10402483
;;; 2.412 to 2.493: 0.12174609
;;; 2.493 to 2.574: 0.1222467
;;; 2.574 to 2.654: 0.13035643
;;; 2.654 to 2.735: 0.13576292
;;; 2.735 to 2.816: 0.09911894
;;; 2.816 to 2.897: 0.06938326
;;; 2.897 to 2.978: 0.044353224
;;; 2.978 to 3.059: 0.0155186225
;;; 3.059 to 3.140: 0.0066079297
;;; 3.140 to 3.221: 0.0

;;; The distribution is now closer to (pi), but not all the way collapsed. Notably, all the probability mass lies
;;; outside the original intervals off of just 1 (forward backward) pair . This is probably due to the  weight
;;; function being tightly tied to the sensor measurement - the standard deviation on a sensor measurement is
;;; roughly 0.314 radians, any state outside that range is likely to be heavily discounted during resampled.
;;;
;;; What if we had a sensor sequence consistent with the intervals [2 pi/3 , 3 pi/4] and [pi/2 pi/3], but an
;;; initial belief centered around pi? That is what our third example focuses on. 


;;;;; **********
;;;;;; OLD MODEL
;;;;; **********
;;;
;;; Left unperturbed so results from parts 1&2 can be checked.
;;;
;;; New model definitions do not overload what's put here.

(defparameter *actions* '(:forward :backward))
(defparameter *states* '(0 1 2 3 4))
(defparameter *sensors* '(:even :odd))

 (defun states () *states*)

;;; A description of our action model


 (defparameter *action-model*
        #3A((  ;; forward
                (.5 .5 0 0 0)
                (0 .5 .5 0 0)
                (0 0 .5 .5 0)
                (0 0 0 .5 .5)
                (.5 0 0 0 .5))

                ;; backward
                ((.5 0 0 0 .5)
                (.5 .5 0 0 0)
                (0 .5 .5 0 0)
                (0 0 .5 .5 0)
                (0 0 0 .5 .5)))
    "A table showing P(new | old, action) where action is either 
forward or backward, old is the row value, and new is the column 
value.  This basically models the notion that, in a toroidal 
one-dimensional world, if you move forward (or backward), you
have a 1/2 probability of accidentally staying put."
        ;;; note that the values in each row must sum to 1
)

(defun action-probability (new-state old-state action)
        "Returns the probability that, given an action and 
old state, the new state will be achieved"
        (aref *action-model* 
                (position action *actions*)
                (position old-state *states*)
                (position new-state *states*)))



;;;; A description of our sensor model


(defparameter *sensor-model*
        #2A(    (0.95 0.05)
                (0.05 0.95)
                (0.95 0.05)
                (0.05 0.95)
                (0.95 0.05))
    "A table showing P(sensor | state), where state is the row 
value, and sensor is the column value.  This basically models 
the notion that in a room of index i (0 to 4), your sensor 3/4 
of the time will give you 0 if i is even and 1 if i is odd.  
I dunno, maybe the rooms are light and dark or something.  :-)"
        ;;; note that the values in each row must sum to 1
)

(defun sensor-probability (sensor state)
        "Returns the probability that, given a state, that the
sensor will be the given value"
        (aref *sensor-model* 
                (position state *states*)
                (position sensor *sensors*)))


;;;;;; Some utility functions for normalizing tables, 
;;;;;; creating particles, and counting
;;;;;; particles at the end.

(defun normalize (lis)
        "Normalizes a table"
        (let ((sum (+ 0.0 (apply #'+ lis)))) ;; force to float
                (mapcar (lambda (elt) (/ elt sum)) lis)))

(defun random-elt (seq)
  "Returns a random element from a sequence"
  (elt seq (random (length seq))))

(defmacro gather (times &rest body)
  "Format:  (gather num-times -body-)

Calls the -body- code num-times times.  Each time the value of the last expression
in -body- is appended to a list.  The full list is then returned."
  (let ((x (gensym)) (y (gensym)))
    `(let (,y) (dotimes (,x ,times) (push (progn ,@body) ,y)) (nreverse ,y))))

(defun counts (distribution)
  "Counts the number of times a state particle appears in the distribution"
  (mapcar (lambda (state) (count state distribution)) *states*))


;;; **********
;;; NEW MODEL
;;; *********
;;;
;;; Continuous/particle filter update
;;;
;;; we need to sample from gaussians, but base lisp only supports uniform random sampling.
;;; looked into the box-muller transform which allows sampling a gaussian with
;;;
;;; std 1 and mean 0
;;;
;;; from two uniform random samples. Just what we need!
(defun box-muller ()
  "Generates (two (we keep one)) random samples from a Gaussian of std 1 centered on 0"
  (let* ((theta (random (* 2.0 pi)))
	 (radius (sqrt (* -2.0 (log (max 1e-15 (random 1.0)))))))
	 (if (= 0 (random 2)) (* radius (cos theta)) (* radius (sin theta)))))

;;; with box-muller, we can sample arbitrary 1-D gaussians. 
(defun sample-from-gaussian (mean std)
  "Gives a sample on a Gaussian (mean, std)"
  (+ mean (* std (box-muller))))

;;; NEW ACTION MODEL
;;;
;;; actions still defined to be forward/backward. We just change what that means. 
;;;
;;; fix a mean "distance" forwards or backwards creates,
;;; as well as a standard error in that motion.
;;;
;;; Note we've chosen the forward motor to be a little stronger -
;;; error in true distance travelled is always proportional to
;;; the mean distance for either direction

(defparameter *2pi* (* 2.0 pi))
(defparameter *forward-mean-angle* (/ *2pi* 8.0))
(defparameter *forward-angle-error* (/ *forward-mean-angle* 5))
(defparameter *backward-mean-angle* (/ (* -1.0 *2pi*) 10.0))
(defparameter *backward-angle-error* (/ *backward-mean-angle* -5.0))

(defparameter *actions-cont* '(:forward :backward))

(defparameter *action-model-cont*
  (list (list *forward-mean-angle* *forward-angle-error*)
   (list *backward-mean-angle* *backward-angle-error*)))

(defun handle-collisions (theta)
  "Handle collisions into the wall during actions"
  
   ;;; here we handle edge cases.
   ;;; for two*pi+k bounce back to 2pi - k
   ;;; for 0-k, bounce back to k
   ;;; edge cases are when a sample is > +4pi or < - 4pi. So long as my
   ;;; standard devaiations for my movement is at most 1/20 * 2pi,
   ;;; I calculated that those samples shouldn't happen during
   ;;; the lifetime of the universe.
   ;;; those edge cases are thus not considerd.
      (cond ((< *2pi* theta) (- (* 2 *2pi*) theta))
	  ((> 0 theta) (* -1.0 theta))
	  (T theta)))

(defun action-result (state action)
  " Given a theta and an action, pick a resulting theta for this action. "
         ;;; this block determines the distribution to sample and samples it.
  (let* ((params (elt *action-model-cont* (position action *actions-cont*)))
	 (motion (first params))
	 (err (second params))
	 (new-guess (+ state motion))
	 (sampled-result (sample-from-gaussian new-guess err)))
    (handle-collisions sampled-result)))
    

;;; NEW SENSOR MODEL
;;;
;;; sensor should give min (theta, 2pi - theta) + epsilon
;;; i.e. distance to wall, with noise
;;; When weighing a particle, it's likelihood is proportional to a gaussian
;;; centered around the sensor reading, i.e. two peaks at theta and 2pi-theta. 

(defparameter *sensor-error* (/ *2pi* 20))

(defun distance-from-wall (state)
  (min state (- *2pi* state)))

(defun sensor-result (state)
  (let ((min-distance (distance-from-wall state)))
    (sample-from-gaussian min-distance *sensor-error*)))

(defun particle-weight (state sensor)
  (exp (- (/ (expt (- (distance-from-wall state) sensor) 2.0) (* 2.0 (expt *sensor-error* 2.0))))))



;;;;; Table-based Bayes Filter.  The function BAYES-FILTER takes a collection of
;;;;; previous beliefs about what state we're in and returns a new collection.
;;;;; Given our simple model above, these beliefs will take the form '(p0 p1 p2 p3 p4),
;;;;; which is just a table of probabilities, one per state.

(defun action-probabilities (new-state action)
  "Given a new state and an action, returns a list of probabilities, one
per old state, of the probability of transitioning to the new state from
that old state given the action"
  (let ((achieve-new-state-probs))
    (loop for old-state in (states)
	  do
	     (push (action-probability new-state old-state action) achieve-new-state-probs))
    (reverse achieve-new-state-probs)))

(defun sensor-probabilities (sensor)
  "Given a sensor value, returns a list of probabilities, one
per state, of the probability of receiving that sensor value given the state"
  (let ((prob-in-state))
    (loop for state in (states)
	  do     
	     (push (sensor-probability sensor state) prob-in-state))
    (reverse prob-in-state)))
     
(defun bayes-filter (action sensor previous-beliefs)
  "Given a previous belief table, an action, and a sensor result, returns a new belief table about what our new states might possibly be. Belief tables are simply lists of the form '(p1 p2 p3 p4 p5 ...) where p1 is the probability for the first state, p2 is the probability for the second state, and so on."
  (let ((belief-table))
    (loop for state in (states)
	  do
	     (push (reduce '+ (mapcar '* previous-beliefs (action-probabilities state action))) belief-table))
    (normalize (mapcar '* (reverse belief-table) (sensor-probabilities sensor)))))

	
;;;;; The particle filter.  The function PARTICLE-FILTER is similar to BAYES-FILTER
;;;;; except that its collections of beliefs take the form of BAGS of STATES.  The same
;;;;; state may repeatedly appear many times in this bag.  The number of times a state is
;;;;; is in the bag basically approximates the probability that we believe we're likely in
;;;;; that state.


;; I use this function to sample a random index from a distribution
;; of the form '(0.1 0.4 0.3 0.2)  for example.  In this example, the
;; returned value would be one of 0, 1, 2, or 3.  You can do this in O(n) if you
;; want, it's okay to be inefficient here.
(defun sample-from-distribution (distribution)
  "Given a possibly non-normalized distribution, in the form of a list of probabilities,
selects from the distribution randomly and returns the index of the selected element."
  (let ((threshold (random (reduce '+ distribution)))
	(sum 0.0)
	(dist-last-index (- (length distribution) 1)))
    (loop for weight in distribution
	  for index from 0 to dist-last-index
	  do
	     (setf sum (+ sum weight))
	     (cond ((>= sum threshold) (return index))
		   ((= index dist-last-index) (return index))))))



;; I use this function to convert a distribution of the form '((x1 w1) (x2 w2) ...)
;; into one like this: '(x1 x1 x1 x2 x2 x3 x4 x4 x4 x5 ...), basically randomly grabbing
;; x's from the previous distribution, proportionally to their weights, and sticking them
;; in the new distribution.  Note that the weights don't have to sum to 1 -- they're
;; weights, not probabilities.
(defun resample-distribution (samples-and-weights &key (sample #'first) (weight #'second))
  "Given a distribution (a list) of M elements, samples M times from that distribution and
forms a new distribution of elements.  Uses a Low-variance 'Stochastic Universal Sampling'
 (wikipedia for it) sampler, which is what the book uses; 
but a plain-old roulette wheel sampler or something else could have
been used just as well.  The function -sample- provides a possible sample in the list, and the function -weight- provides the probability that it should be selected.  By default these are #'first and #'second, presuming that the provided distribution is
of the form ((sample1 weight1) (sample2 weight2) ...)."
  ; following the wikipedia pseudocode
  (let* ((total-fitness (reduce '+ (mapcar weight samples-and-weights)))
	 (num-samples (length samples-and-weights))
	 (distance-between-pointers (/ total-fitness num-samples))
	 (start (random distance-between-pointers))
	 (pointers (loop for i from 0 to (- num-samples 1)
			 collect (+ start (* i distance-between-pointers))))
	 (bag)
	 (fitness-sum 0.0)
	 (point-index 0))
    ; only summing throught the weights once, instead of a new time for each particle. 
    (loop for s-and-w in samples-and-weights
    	  do
	     (setf fitness-sum (+ fitness-sum (funcall weight s-and-w)))
	     (loop while (and (< point-index num-samples) (>= fitness-sum (elt pointers point-index)))
		   do
		      (push (funcall sample s-and-w) bag)
		      (setf point-index (+ point-index 1))))
    bag))
	 

;; *OLD MODEL*
;; Here I grab a random new state selected from the distribution of old ones.
;; note that I'm grabbing from a table -- that doesn't have to be the case,
;; I'm doing it to be compatable with my bayes filter examples.  In fact, you
;; don't have to have a "distribution" at all -- you can just create a function
;; which picks a random new state given the old state and action.  Often these
;; functions are much easier to write than generating a whole probability distribution.
;; And they allow arbitrarily complex continuous distributions to boot.  But here
;; we're going with the simple backward-compatable code...
(defun select-from-action-probabilities (old-state action)
  "Given an old-state and an action, selects a new-state at random and returns it
given the probability distribution P(new-state | old-state, action)."
  (let ((probs)
	(states (states)))
    (loop for new-state in states
	  do
	     (push (action-probability new-state old-state action) probs))
    (elt states (sample-from-distribution (reverse probs)))))



;;;
;;;
;;; *************************
;;; TOGGLE BETWEEN MODELS HERE
;;; *************************
;;;
;;; The old model has discretized states/sensor readings.
;;; That does not play well with the new continuous regime.
;;; We can abstract away *all* the changes to these two wrapper
;;; functions though, which is nice.
;;;
;;; To run the discrete examples, uncomment the OLD MODEL lines
;;;
;;; To try out the continuous examples, uncomment the NEW MODEL lines.
;;;
(defun action-transition (old-state action)
  "given some theta and an action, get a new particle"
  ;;; *NEW MODEL*
  (action-result old-state action))
  ;;; *OLD MODEL*
  ;;; (select-from-action-probabilities old-state action))

(defun weight-function (state sensor)
  "given some particle and a sensor reading, determine its fitness."

  ;;; *NEW MODEL*
  (list state (particle-weight state sensor)))
  ;;; *OLD MODEL*
  ;;; (list state (elt (sensor-probabilities sensor) state)))

;; Now we just do the belief update.  Note that beliefs now are different than they
;; were in the past: they're particles, each one representing a state.  Also note that
;; just as you could have sritten select-from-action-probabilities above without even
;; *having* a distribution, you can do the same for sensor-probability.  That's a major
;; strength of the particle filter.  It makes for easier representations of your
;; "distributions".
(defun particle-filter (action sensor previous-beliefs)
  "Given a previous belief table, an action, and a sensor 
result, returns a  new belief table about what our new states 
might possibly be.  Belief tables are lists of particles.  Each
particle is simply a state."
  (let* ((particles-from-action (mapcar #'(lambda (p) (action-transition p action)) previous-beliefs))
	 (particle-weight-pairs (mapcar #'(lambda (p) (weight-function p sensor)) particles-from-action)))
    (resample-distribution particle-weight-pairs)))

    
  
;;;; SOME TEST FUNCTIONS
;;;; (example-1-bayes)  and (example-1-particle) should return approximately
;;;; the same numbers.  Likewise for (example-2-bayes) and (example-2-particle)

;;;; For convenience:
(defvar *b*)

(defun counts-general (distribution bin-numbers)
  "Expanding on the `counts` function to work on any list of bin-numbers."
  (mapcar (lambda (sample) (count sample distribution)) bin-numbers))  

(defun make-hist (distribution min max num-bins)
  "Returns a normalized probability distribution binned to the caller's specifications."
  (let* ((bin-width (/ (- max min) num-bins))
	 (bin-numbers (loop for i from 0 to (- num-bins 1) collect i))
	 (bins (mapcar #'(lambda (number) (+ min (* number bin-width))) bin-numbers))
	 (bag))
    (loop for sample in distribution
	  do
	     (push (floor (/ (- sample min) bin-width)) bag))
    (list (normalize (counts-general bag bin-numbers)) (nconc bins (list max)))))

(defun print-hist (hist bins)
  "Prints a histogram to the screen in a readable format."
  (format t "~%")
  (format t "~%")
  (format t "Probability mass on [~,3F,~,3F]~%" (first bins) (car (last bins)))
  (format t "~%")
  (loop for count in hist
	for index from 0 to (- (length bins) 2)
	do
	   (format t "~,3F to ~,3F: ~F~%" (elt bins index) (elt bins (+ 1 index)) count))
  (format t "~%")
  (format t "~%"))

(defun generate-particles-on-interval (theta-min theta-max num-particles)
  "Generates a uniform distributions of particles on a provided angle interval"
  (let ((bag))
    (dotimes (iter num-particles)
      (push (+ theta-min (random (- theta-max theta-min))) bag))
    bag))


;;; NEW TEST FUNCTIONS

(defun sensor-vals-example-1 (action-sequence)
  "Give average sensor readings for the action sequence outlined in example-1-continuous"
  (let ((bag)
	(current-theta pi))
    (loop for action in action-sequence
	  do
	     (if (eq action :forward) (setf current-theta (+ current-theta *forward-mean-angle*))
		 (setf current-theta (+ current-theta *backward-mean-angle*)))
	     (push (distance-from-wall current-theta) bag))
    (reverse bag)))


(defun example-1-continuous ()
  "First example calculation, simulates a robot that starts opposite the wall, near theta=pi"
  (let* ((particles (generate-particles-on-interval (- pi 0.2) (+ pi 0.2) 10000))
	(action-sequence (list :forward :backward :forward :backward 
			   :forward :backward :forward :backward
			   :forward :backward))
	(sensor-sequence (sensor-vals-example-1 action-sequence)))
    (loop for action in action-sequence
	  for sensor in sensor-sequence
	  do
	     (setf particles (particle-filter action sensor particles)))
    (apply #'print-hist (make-hist particles (reduce #'min particles) (reduce #'max particles) 20.0))))

(defun example-2-continuous ()
  "Second example calculation, simulates a robot with a mistaken initial belief."
  (let* ((particles (append
		     (generate-particles-on-interval (/ (* 2 pi) 3) (/ (* 3 pi) 4) 5000)
		     (generate-particles-on-interval (/ pi 4) (/ pi 3) 5000)))
	 ;;;
	;;; SAME action/senor sequence, much different initial belief	    
;;;	 (action-sequence (list :forward :backward :forward :backward 
;;;			       :forward :backward :forward :backward
;;;			       :forward :backward))
	 (action-sequence (list :forward))
	 (sensor-sequence (sensor-vals-example-1 action-sequence)))
    (loop for action in action-sequence
	  for sensor in sensor-sequence
	  do
	     (setf particles (particle-filter action sensor particles)))
    (apply #'print-hist (make-hist particles (reduce #'min particles) (reduce #'max particles) 20.0))))
    
;;;; Example 1 : move forward through the environment,
;;;; then back to where you started.  You don't know where
;;;; you started.  Where are you likely to be?


;;; with the bayes filter
(defun example-1-bayes ()
  (let ((b (normalize '(1 1 1 1 1))))
    (setf b (bayes-filter :forward :odd b))
    (setf b (bayes-filter :forward :even b)) 
    (setf b (bayes-filter :forward :odd b)) 
    (setf b (bayes-filter :forward :even b)) 
    (setf b (bayes-filter :forward :even b)) 
    (setf b (bayes-filter :forward :odd b))
    (setf b (bayes-filter :forward :even b))
    (setf b (bayes-filter :forward :odd b))
    (setf b (bayes-filter :backward :even b)) 
    (setf b (bayes-filter :backward :odd b)) 
    (setf b (bayes-filter :backward :even b)) 
    (setf b (bayes-filter :backward :even b)) 
    (setf b (bayes-filter :backward :odd b)) 
    (setf b (bayes-filter :backward :even b)) 
    (setf b (bayes-filter :backward :odd b))
    (setf b (bayes-filter :backward :even b))
    (format t "Bayes Filter Results: ~a" b)))


;;; with the particle filter
(defun example-1-particle ()
  (let ((b (gather 10000 (random-elt *states*))))
    (setf b (particle-filter :forward :odd b))
    (print "1/16")
    (setf b (particle-filter :forward :even b))
    (print "2/16")
    (setf b (particle-filter :forward :odd b))
    (print "3/16")
    (setf b (particle-filter :forward :even b))
    (print "4/16")
    (setf b (particle-filter :forward :even b))
    (print "5/16")
    (setf b (particle-filter :forward :odd b))
    (print "6/16")
    (setf b (particle-filter :forward :even b))
    (print "7/16")
    (setf b (particle-filter :forward :odd b))
    (print "8/16")
    (setf b (particle-filter :backward :even b))
    (print "9/16")
    (setf b (particle-filter :backward :odd b))
    (print "10/16")
    (setf b (particle-filter :backward :even b))
    (print "11/16")
    (setf b (particle-filter :backward :even b))
    (print "12/16")
    (setf b (particle-filter :backward :odd b))
    (print "13/16")
    (setf b (particle-filter :backward :even b))
    (print "14/16")
    (setf b (particle-filter :backward :odd b))
    (print "15/16")
    (setf b (particle-filter :backward :even b))
    (print "16/16")
    (format t "Particle Filter Results: ~a" (normalize (counts b)))))





;;; example two: If you knew you were in room 0,
;;; then moved forward twice and backward thrice,
;;; and your sensors said even odd even even odd,
;;; where you you likely to be?

;;; with the bayes filter

(defun example-2-bayes ()
  (let ((b '(1 0 0 0 0)))
    (setf b (bayes-filter :forward :even b))
    (setf b (bayes-filter :forward :odd b))
    (setf b (bayes-filter :backward :even b))
    (setf b (bayes-filter :backward :even b))
    (setf b (bayes-filter :backward :odd b))
    (format t "Bayes Filter Results: ~a" b)))



;;; with the particle filter
(defun example-2-particle ()
  (let ((b (gather 10000 0)))
    (setf b (particle-filter :forward :even b))
    (setf b (particle-filter :forward :odd b))
    (setf b (particle-filter :backward :even b))
    (setf b (particle-filter :backward :even b))
    (setf b (particle-filter :backward :odd b))
    (format t "Particle Filter Results: ~a" (normalize (counts b)))))
