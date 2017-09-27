#lang racket

(require "declarations.rkt")
(provide buildTree calcForces moveparticles)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (buildTree initialArea particles)
  (cond [(= (length particles) 1) (gnode (particle-mass (car particles)) (particle-posn (car particles)) '())]
        [(= (length particles) 0) (gnode 0 (sq-cen initialArea) '())]
        [else (let* ([1st-quad-area (1st-ar initialArea)]
                     [2nd-quad-area (2nd-ar initialArea)]
                     [3rd-quad-area (3rd-ar initialArea)]
                     [4th-quad-area (4th-ar initialArea)])
                (gnode (particle-mass (cm particles)) (particle-posn (cm particles))
                       (list (buildTree 2nd-quad-area (quad-par 2nd-quad-area particles))
                             (buildTree 1st-quad-area (quad-par 1st-quad-area particles))
                             (buildTree 3rd-quad-area (quad-par 3rd-quad-area particles))
                             (buildTree 4th-quad-area (quad-par 4th-quad-area particles)))))]))

(define (cm particles)
  (let* [(mass-list (lc (particle-mass p) : p <- particles))
         (m (foldr + 0 mass-list))
         (pos-x (f mass-list (lc (vec-x (particle-posn p)) : p <- particles)))
         (pos-y (f mass-list (lc (vec-y (particle-posn p)) : p <- particles)))
         (vel-x (f mass-list (lc (vec-x (particle-velocity p)) : p <- particles)))
         (vel-y (f mass-list (lc (vec-y (particle-velocity p)) : p <- particles)))]
  (particle m (vec pos-x pos-y) (vec vel-x vel-y))))

(define (f l1 l2)
  (if (or (null? l1) (null? l2)) 0
      (let* ([m1 (car l1)]
             [m2 (foldr + 0 (cdr l1))])
        (/ (+ (* m1 (car l2)) (* m2 (f (cdr l1) (cdr l2)))) (+ m1 m2)))))

(define (sq-cen box)
  (vec (/ (+ (bbox-llx box) (bbox-rux box)) 2) (/ (+ (bbox-lly box) (bbox-ruy box)) 2)))

(define (1st-ar initialArea)
  (bbox (/ (+ (bbox-llx initialArea) (bbox-rux initialArea)) 2)
        (/ (+ (bbox-lly initialArea) (bbox-ruy initialArea)) 2)
        (bbox-rux initialArea)
        (bbox-ruy initialArea)))

(define (2nd-ar initialArea)
  (bbox (bbox-llx initialArea)
        (/ (+ (bbox-lly initialArea) (bbox-ruy initialArea)) 2)
        (/ (+ (bbox-llx initialArea) (bbox-rux initialArea)) 2)
        (bbox-ruy initialArea)))

(define (3rd-ar initialArea)
  (bbox (bbox-llx initialArea)
        (bbox-lly initialArea)
        (/ (+ (bbox-llx initialArea) (bbox-rux initialArea)) 2)
        (/ (+ (bbox-lly initialArea) (bbox-ruy initialArea)) 2)))

(define (4th-ar initialArea)
  (bbox (/ (+ (bbox-llx initialArea) (bbox-rux initialArea)) 2)
        (bbox-lly initialArea)
        (bbox-rux initialArea)
        (/ (+ (bbox-lly initialArea) (bbox-ruy initialArea)) 2)))

(define (quad-par box particles)
  (if (null? particles) '()
      (let* ([par (car particles)]
             [posx (vec-x (particle-posn par))]
             [posy (vec-y (particle-posn par))])
        (if (and (> posx (bbox-llx box)) (< posx (bbox-rux box))
                 (> posy (bbox-lly box)) (< posy (bbox-ruy box))) (cons par (quad-par box (cdr particles)))
                                                                  (quad-par box (cdr particles))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (calcForces initialArea tree particles)
  (if (null? particles) '()
      (cons (force initialArea (car particles) tree) (calcForces initialArea tree (cdr particles)))))

(define (force initialArea par tree)
  (if (= (particle-mass par) 0) (vec 0 0)
      (let* ([m1 (particle-mass par)]
             [pos1 (particle-posn par)]
             [m2 (gnode-mass tree)]
             [pos2 (gnode-posn tree)]
             [new-area (1st-ar initialArea)])
        (if (or (> (/ (distance pos1 pos2) (- (bbox-rux initialArea) (bbox-llx initialArea))) theta)
                (null? (gnode-subtrees tree))) (force1 m1 pos1 m2 pos2)
            (vec-sum (force new-area par (car (gnode-subtrees tree)))
                     (force new-area par (cadr (gnode-subtrees tree)))
                     (force new-area par (caddr (gnode-subtrees tree)))
                     (force new-area par (cadddr (gnode-subtrees tree))))))))

(define (distance pos1 pos2)
  (sqrt (+ (* (- (vec-x pos1) (vec-x pos2)) (- (vec-x pos1) (vec-x pos2)))
           (* (- (vec-y pos1) (vec-y pos2)) (- (vec-y pos1) (vec-y pos2))))))

(define (force1 m1 pos1 m2 pos2)
  (if (= 0 (distance pos1 pos2)) (vec 0 0)
      (let* ([d (distance pos1 pos2)]
             [fx (/ (* m1 m2 g (- (vec-x pos2) (vec-x pos1))) (* d d d))]
             [fy (/ (* m1 m2 g (- (vec-y pos2) (vec-y pos1))) (* d d d))])
        (vec fx fy))))

(define (vec-sum a b c d)
  (vec (+ (vec-x a) (vec-x b) (vec-x c) (vec-x d)) (+ (vec-y a) (vec-y b) (vec-y c) (vec-y d))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (moveparticles particles forces)
  (if (null? particles) '()
      (cons (move-particle (car particles) (car forces)) (moveparticles (cdr particles) (cdr forces)))))

(define (move-particle par forc)
  (let* ([t timeslice]
         [a (vec (/ (vec-x forc) (particle-mass par)) (/ (vec-y forc) (particle-mass par)))])
    (new-particle par a t)))

(define (new-particle par a t)
  (let* ([new-px (+ (vec-x (particle-posn par)) (* (vec-x (particle-velocity par)) t) (* 0.5 (vec-x a) t t))]
         [new-py (+ (vec-y (particle-posn par)) (* (vec-y (particle-velocity par)) t) (* 0.5 (vec-y a) t t))]
         [new-vx (+ (vec-x (particle-velocity par)) (* (vec-x a) t))]
         [new-vy (+ (vec-y (particle-velocity par)) (* (vec-y a) t))])
    (particle (particle-mass par) (vec new-px new-py) (vec new-vx new-vy))))