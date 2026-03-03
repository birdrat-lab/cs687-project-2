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
  (let* ((particles-from-action (mapcar #'(lambda (p) (select-from-action-probabilities p action)) previous-beliefs))
	 (sensor-probs (sensor-probabilities sensor))
	 (particle-weight-pairs (mapcar #'(lambda (p) (list p (elt sensor-probs p))) particles-from-action)))
    (resample-distribution particle-weight-pairs)))


;;;; SOME TEST FUNCTIONS
;;;; (example-1-bayes)  and (example-1-particle) should return approximately
;;;; the same numbers.  Likewise for (example-2-bayes) and (example-2-particle)

;;;; For convenience:
(defvar *b*)


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
    (setf b (particle-filter :forward :even b)) 
    (setf b (particle-filter :forward :odd b)) 
    (setf b (particle-filter :forward :even b)) 
    (setf b (particle-filter :forward :even b)) 
    (setf b (particle-filter :forward :odd b))
    (setf b (particle-filter :forward :even b))
    (setf b (particle-filter :forward :odd b))
    (setf b (particle-filter :backward :even b)) 
    (setf b (particle-filter :backward :odd b)) 
    (setf b (particle-filter :backward :even b)) 
    (setf b (particle-filter :backward :even b)) 
    (setf b (particle-filter :backward :odd b)) 
    (setf b (particle-filter :backward :even b)) 
    (setf b (particle-filter :backward :odd b))
    (setf b (particle-filter :backward :even b))
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
