% Calculate inertia matrix for parallel robot
% P3RRRRR2G3P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR2G3P3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:11:57
% EndTime: 2020-03-09 21:11:58
% DurationCPUTime: 1.03s
% Computational Cost: add. (1743->198), mult. (3279->407), div. (792->11), fcn. (3042->27), ass. (0->217)
t494 = sin(qJ(2,3));
t616 = pkin(1) * t494;
t497 = sin(qJ(2,2));
t615 = pkin(1) * t497;
t500 = sin(qJ(2,1));
t614 = pkin(1) * t500;
t503 = cos(qJ(2,3));
t613 = pkin(1) * t503;
t504 = cos(qJ(1,3));
t612 = pkin(1) * t504;
t506 = cos(qJ(2,2));
t611 = pkin(1) * t506;
t507 = cos(qJ(1,2));
t610 = pkin(1) * t507;
t509 = cos(qJ(2,1));
t609 = pkin(1) * t509;
t510 = cos(qJ(1,1));
t608 = pkin(1) * t510;
t607 = Ifges(3,1) + Ifges(2,3);
t493 = sin(qJ(3,3));
t606 = Ifges(3,4) * t493;
t496 = sin(qJ(3,2));
t605 = Ifges(3,4) * t496;
t499 = sin(qJ(3,1));
t604 = Ifges(3,4) * t499;
t490 = legFrame(3,2);
t469 = sin(t490);
t472 = cos(t490);
t581 = t472 * t493;
t495 = sin(qJ(1,3));
t448 = t495 * t494 - t504 * t503;
t502 = cos(qJ(3,3));
t594 = t448 * t502;
t433 = t469 * t594 + t581;
t480 = 0.1e1 / t502;
t603 = t433 * t480;
t587 = t469 * t493;
t434 = -t472 * t594 + t587;
t602 = t434 * t480;
t491 = legFrame(2,2);
t470 = sin(t491);
t473 = cos(t491);
t579 = t473 * t496;
t498 = sin(qJ(1,2));
t449 = t498 * t497 - t507 * t506;
t505 = cos(qJ(3,2));
t593 = t449 * t505;
t435 = t470 * t593 + t579;
t483 = 0.1e1 / t505;
t601 = t435 * t483;
t585 = t470 * t496;
t436 = -t473 * t593 + t585;
t600 = t436 * t483;
t492 = legFrame(1,2);
t471 = sin(t492);
t474 = cos(t492);
t577 = t474 * t499;
t501 = sin(qJ(1,1));
t450 = t501 * t500 - t510 * t509;
t508 = cos(qJ(3,1));
t592 = t450 * t508;
t437 = t471 * t592 + t577;
t486 = 0.1e1 / t508;
t599 = t437 * t486;
t583 = t471 * t499;
t438 = -t474 * t592 + t583;
t598 = t438 * t486;
t442 = t502 * (-mrSges(3,2) * t616 + Ifges(3,6)) - t493 * (mrSges(3,1) * t616 - Ifges(3,5));
t597 = t442 * t480;
t443 = t505 * (-mrSges(3,2) * t615 + Ifges(3,6)) - t496 * (mrSges(3,1) * t615 - Ifges(3,5));
t596 = t443 * t483;
t444 = t508 * (-mrSges(3,2) * t614 + Ifges(3,6)) - t499 * (mrSges(3,1) * t614 - Ifges(3,5));
t595 = t444 * t486;
t466 = sin(qJ(1,3) + qJ(2,3));
t476 = 0.1e1 / t494;
t591 = t466 * t476;
t467 = sin(qJ(1,2) + qJ(2,2));
t477 = 0.1e1 / t497;
t590 = t467 * t477;
t468 = sin(qJ(1,1) + qJ(2,1));
t478 = 0.1e1 / t500;
t589 = t468 * t478;
t588 = t469 * t480;
t586 = t470 * t483;
t584 = t471 * t486;
t582 = t472 * t480;
t580 = t473 * t483;
t578 = t474 * t486;
t576 = t476 * t480;
t481 = 0.1e1 / t502 ^ 2;
t575 = t476 * t481;
t512 = 0.1e1 / pkin(1);
t574 = t476 * t512;
t573 = t477 * t483;
t484 = 0.1e1 / t505 ^ 2;
t572 = t477 * t484;
t571 = t477 * t512;
t570 = t478 * t486;
t487 = 0.1e1 / t508 ^ 2;
t569 = t478 * t487;
t568 = t478 * t512;
t511 = 0.1e1 / pkin(2);
t567 = t480 * t511;
t566 = t481 * t511;
t565 = t483 * t511;
t564 = t484 * t511;
t563 = t486 * t511;
t562 = t487 * t511;
t561 = 0.2e1 * t606;
t560 = 0.2e1 * t605;
t559 = 0.2e1 * t604;
t558 = t493 * t613;
t557 = t496 * t611;
t556 = t499 * t609;
t479 = t502 ^ 2;
t555 = pkin(2) * t448 * t479;
t482 = t505 ^ 2;
t554 = pkin(2) * t449 * t482;
t485 = t508 ^ 2;
t553 = pkin(2) * t450 * t485;
t488 = Ifges(3,2) - Ifges(3,1);
t457 = t488 * t479;
t552 = t457 + t607;
t458 = t488 * t482;
t551 = t458 + t607;
t459 = t488 * t485;
t550 = t459 + t607;
t421 = -t469 * t555 + (-pkin(2) * t581 + t469 * t612) * t502 - t472 * t558;
t549 = t421 * t575;
t422 = -t470 * t554 + (-pkin(2) * t579 + t470 * t610) * t505 - t473 * t557;
t548 = t422 * t572;
t423 = -t471 * t553 + (-pkin(2) * t577 + t471 * t608) * t508 - t474 * t556;
t547 = t423 * t569;
t424 = t472 * t555 + (-pkin(2) * t587 - t472 * t612) * t502 - t469 * t558;
t546 = t424 * t575;
t425 = t473 * t554 + (-pkin(2) * t585 - t473 * t610) * t505 - t470 * t557;
t545 = t425 * t572;
t426 = t474 * t553 + (-pkin(2) * t583 - t474 * t608) * t508 - t471 * t556;
t544 = t426 * t569;
t463 = mrSges(3,1) * t613;
t489 = mrSges(2,2) - mrSges(3,3);
t515 = (-(t493 * mrSges(3,2) - mrSges(2,1)) * t503 - t489 * t494) * pkin(1);
t430 = (t463 + t561) * t502 + t515 + t552;
t543 = t430 * t566;
t464 = mrSges(3,1) * t611;
t514 = (-(t496 * mrSges(3,2) - mrSges(2,1)) * t506 - t489 * t497) * pkin(1);
t431 = (t464 + t560) * t505 + t514 + t551;
t542 = t431 * t564;
t465 = mrSges(3,1) * t609;
t513 = (-(t499 * mrSges(3,2) - mrSges(2,1)) * t509 - t489 * t500) * pkin(1);
t432 = (t465 + t559) * t508 + t513 + t550;
t541 = t432 * t562;
t540 = t433 * t576;
t539 = t434 * t576;
t538 = t435 * t573;
t537 = t436 * t573;
t536 = t437 * t570;
t535 = t438 * t570;
t439 = pkin(2) * (t504 * t494 + t495 * t503) * t502 + t495 * pkin(1);
t534 = t439 * t576;
t533 = t439 * t567;
t440 = pkin(2) * (t507 * t497 + t498 * t506) * t505 + t498 * pkin(1);
t532 = t440 * t573;
t531 = t440 * t565;
t441 = pkin(2) * (t510 * t500 + t501 * t509) * t508 + t501 * pkin(1);
t530 = t441 * t570;
t529 = t441 * t563;
t445 = t502 * t561 + t552;
t528 = t445 * t566;
t446 = t505 * t560 + t551;
t527 = t446 * t564;
t447 = t508 * t559 + t550;
t526 = t447 * t562;
t451 = Ifges(3,5) * t493 + Ifges(3,6) * t502;
t525 = t451 * t566;
t452 = Ifges(3,5) * t496 + Ifges(3,6) * t505;
t524 = t452 * t564;
t453 = Ifges(3,5) * t499 + Ifges(3,6) * t508;
t523 = t453 * t562;
t522 = t469 * t567;
t521 = t470 * t565;
t520 = t471 * t563;
t519 = t472 * t567;
t518 = t473 * t565;
t517 = t474 * t563;
t516 = Ifges(1,3) + (m(2) + m(3)) * pkin(1) ^ 2 + t607;
t429 = t459 + 0.2e1 * (t465 + t604) * t508 + 0.2e1 * t513 + t516;
t428 = t458 + 0.2e1 * (t464 + t605) * t505 + 0.2e1 * t514 + t516;
t427 = t457 + 0.2e1 * (t463 + t606) * t502 + 0.2e1 * t515 + t516;
t420 = (-t444 * t468 + t453 * t529) * t568;
t419 = (-t443 * t467 + t452 * t531) * t571;
t418 = (-t442 * t466 + t451 * t533) * t574;
t417 = (-t432 * t468 + t447 * t529) * t568;
t416 = (-t431 * t467 + t446 * t531) * t571;
t415 = (-t430 * t466 + t445 * t533) * t574;
t414 = (-t429 * t468 + t432 * t529) * t568;
t413 = (-t428 * t467 + t431 * t531) * t571;
t412 = (-t427 * t466 + t430 * t533) * t574;
t411 = Ifges(3,3) * t520 + (t426 * t523 + t438 * t595) * t568;
t410 = Ifges(3,3) * t517 + (t423 * t523 + t437 * t595) * t568;
t409 = Ifges(3,3) * t521 + (t425 * t524 + t436 * t596) * t571;
t408 = Ifges(3,3) * t518 + (t422 * t524 + t435 * t596) * t571;
t407 = Ifges(3,3) * t522 + (t424 * t525 + t434 * t597) * t574;
t406 = Ifges(3,3) * t519 + (t421 * t525 + t433 * t597) * t574;
t405 = t453 * t520 + (t426 * t526 + t432 * t598) * t568;
t404 = t453 * t517 + (t423 * t526 + t432 * t599) * t568;
t403 = t452 * t521 + (t425 * t527 + t431 * t600) * t571;
t402 = t452 * t518 + (t422 * t527 + t431 * t601) * t571;
t401 = t451 * t522 + (t424 * t528 + t430 * t602) * t574;
t400 = t451 * t519 + (t421 * t528 + t430 * t603) * t574;
t399 = t444 * t520 + (t426 * t541 + t429 * t598) * t568;
t398 = t444 * t517 + (t423 * t541 + t429 * t599) * t568;
t397 = t443 * t521 + (t425 * t542 + t428 * t600) * t571;
t396 = t443 * t518 + (t422 * t542 + t428 * t601) * t571;
t395 = t442 * t522 + (t424 * t543 + t427 * t602) * t574;
t394 = t442 * t519 + (t421 * t543 + t427 * t603) * t574;
t1 = [m(4) + (t395 * t539 + t397 * t537 + t399 * t535) * t512 + (t407 * t588 + t409 * t586 + t411 * t584 + (t401 * t546 + t403 * t545 + t405 * t544) * t512) * t511, (t395 * t540 + t397 * t538 + t399 * t536) * t512 + (t407 * t582 + t409 * t580 + t411 * t578 + (t401 * t549 + t403 * t548 + t405 * t547) * t512) * t511, (-t395 * t591 - t397 * t590 - t399 * t589 + (t401 * t534 + t403 * t532 + t405 * t530) * t511) * t512; (t394 * t539 + t396 * t537 + t398 * t535) * t512 + (t406 * t588 + t408 * t586 + t410 * t584 + (t400 * t546 + t402 * t545 + t404 * t544) * t512) * t511, m(4) + (t394 * t540 + t396 * t538 + t398 * t536) * t512 + (t406 * t582 + t408 * t580 + t410 * t578 + (t400 * t549 + t402 * t548 + t404 * t547) * t512) * t511, (-t394 * t591 - t396 * t590 - t398 * t589 + (t400 * t534 + t402 * t532 + t404 * t530) * t511) * t512; (t412 * t539 + t413 * t537 + t414 * t535) * t512 + (t418 * t588 + t419 * t586 + t420 * t584 + (t415 * t546 + t416 * t545 + t417 * t544) * t512) * t511, (t412 * t540 + t413 * t538 + t414 * t536) * t512 + (t418 * t582 + t419 * t580 + t420 * t578 + (t415 * t549 + t416 * t548 + t417 * t547) * t512) * t511, m(4) + (-t412 * t591 - t413 * t590 - t414 * t589 + (t415 * t534 + t416 * t532 + t417 * t530) * t511) * t512;];
MX  = t1;
