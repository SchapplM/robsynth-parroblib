% Calculate inertia matrix for parallel robot
% P3RRRRR2G2A0
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
% Datum: 2020-03-09 21:10
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR2G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:08:49
% EndTime: 2020-03-09 21:08:50
% DurationCPUTime: 1.00s
% Computational Cost: add. (1974->213), mult. (3246->396), div. (759->14), fcn. (2943->39), ass. (0->218)
t641 = -2 * pkin(1);
t512 = sin(qJ(2,3));
t640 = pkin(1) * t512;
t513 = sin(qJ(1,3));
t639 = pkin(1) * t513;
t515 = sin(qJ(2,2));
t638 = pkin(1) * t515;
t516 = sin(qJ(1,2));
t637 = pkin(1) * t516;
t518 = sin(qJ(2,1));
t636 = pkin(1) * t518;
t519 = sin(qJ(1,1));
t635 = pkin(1) * t519;
t521 = cos(qJ(2,3));
t634 = pkin(1) * t521;
t524 = cos(qJ(2,2));
t633 = pkin(1) * t524;
t527 = cos(qJ(2,1));
t632 = pkin(1) * t527;
t631 = Ifges(3,1) + Ifges(2,3);
t511 = sin(qJ(3,3));
t630 = Ifges(3,4) * t511;
t514 = sin(qJ(3,2));
t629 = Ifges(3,4) * t514;
t517 = sin(qJ(3,1));
t628 = Ifges(3,4) * t517;
t508 = legFrame(3,2);
t487 = sin(t508);
t490 = cos(t508);
t602 = t490 * t511;
t522 = cos(qJ(1,3));
t466 = t522 * t512 + t513 * t521;
t520 = cos(qJ(3,3));
t615 = t466 * t520;
t451 = -t487 * t615 + t602;
t498 = 0.1e1 / t520;
t627 = t451 * t498;
t608 = t487 * t511;
t452 = t490 * t615 + t608;
t626 = t452 * t498;
t509 = legFrame(2,2);
t488 = sin(t509);
t491 = cos(t509);
t600 = t491 * t514;
t525 = cos(qJ(1,2));
t467 = t525 * t515 + t516 * t524;
t523 = cos(qJ(3,2));
t614 = t467 * t523;
t453 = -t488 * t614 + t600;
t501 = 0.1e1 / t523;
t625 = t453 * t501;
t606 = t488 * t514;
t454 = t491 * t614 + t606;
t624 = t454 * t501;
t510 = legFrame(1,2);
t489 = sin(t510);
t492 = cos(t510);
t598 = t492 * t517;
t528 = cos(qJ(1,1));
t468 = t528 * t518 + t519 * t527;
t526 = cos(qJ(3,1));
t613 = t468 * t526;
t455 = -t489 * t613 + t598;
t504 = 0.1e1 / t526;
t623 = t455 * t504;
t604 = t489 * t517;
t456 = t492 * t613 + t604;
t622 = t456 * t504;
t457 = (-mrSges(3,2) * t640 + Ifges(3,6)) * t520 - t511 * (mrSges(3,1) * t640 - Ifges(3,5));
t621 = t457 * t498;
t458 = (-mrSges(3,2) * t638 + Ifges(3,6)) * t523 - t514 * (mrSges(3,1) * t638 - Ifges(3,5));
t620 = t458 * t501;
t459 = (-mrSges(3,2) * t636 + Ifges(3,6)) * t526 - t517 * (mrSges(3,1) * t636 - Ifges(3,5));
t619 = t459 * t504;
t577 = qJ(2,3) + qJ(3,3);
t578 = qJ(2,3) - qJ(3,3);
t618 = (t522 * t641 + (-cos(qJ(1,3) + t578) - cos(qJ(1,3) + t577)) * pkin(2)) / (sin(t577) + sin(t578));
t579 = qJ(2,2) + qJ(3,2);
t580 = qJ(2,2) - qJ(3,2);
t617 = (t525 * t641 + (-cos(qJ(1,2) + t580) - cos(qJ(1,2) + t579)) * pkin(2)) / (sin(t579) + sin(t580));
t581 = qJ(2,1) + qJ(3,1);
t582 = qJ(2,1) - qJ(3,1);
t616 = (t528 * t641 + (-cos(qJ(1,1) + t582) - cos(qJ(1,1) + t581)) * pkin(2)) / (sin(t581) + sin(t582));
t494 = 0.1e1 / t512;
t612 = cos(qJ(1,3) + qJ(2,3)) * t494;
t495 = 0.1e1 / t515;
t611 = cos(qJ(1,2) + qJ(2,2)) * t495;
t496 = 0.1e1 / t518;
t610 = cos(qJ(1,1) + qJ(2,1)) * t496;
t609 = t487 * t498;
t607 = t488 * t501;
t605 = t489 * t504;
t603 = t490 * t498;
t601 = t491 * t501;
t599 = t492 * t504;
t597 = t494 * t498;
t499 = 0.1e1 / t520 ^ 2;
t596 = t494 * t499;
t530 = 1 / pkin(1);
t595 = t494 * t530;
t594 = t495 * t501;
t502 = 0.1e1 / t523 ^ 2;
t593 = t495 * t502;
t592 = t495 * t530;
t591 = t496 * t504;
t505 = 0.1e1 / t526 ^ 2;
t590 = t496 * t505;
t589 = t496 * t530;
t529 = 0.1e1 / pkin(2);
t588 = t498 * t529;
t587 = t499 * t529;
t586 = t501 * t529;
t585 = t502 * t529;
t584 = t504 * t529;
t583 = t505 * t529;
t576 = 0.2e1 * t630;
t575 = 0.2e1 * t629;
t574 = 0.2e1 * t628;
t497 = t520 ^ 2;
t573 = t466 * t497 * pkin(2);
t500 = t523 ^ 2;
t572 = t467 * t500 * pkin(2);
t503 = t526 ^ 2;
t571 = t468 * t503 * pkin(2);
t570 = t511 * t634;
t569 = t514 * t633;
t568 = t517 * t632;
t506 = Ifges(3,2) - Ifges(3,1);
t478 = t506 * t497;
t567 = t478 + t631;
t479 = t506 * t500;
t566 = t479 + t631;
t480 = t506 * t503;
t565 = t480 + t631;
t439 = t487 * t573 + (-pkin(2) * t602 + t487 * t639) * t520 - t490 * t570;
t564 = t439 * t596;
t440 = -t490 * t573 + (-pkin(2) * t608 - t490 * t639) * t520 - t487 * t570;
t563 = t440 * t596;
t441 = t488 * t572 + (-pkin(2) * t600 + t488 * t637) * t523 - t491 * t569;
t562 = t441 * t593;
t442 = -t491 * t572 + (-pkin(2) * t606 - t491 * t637) * t523 - t488 * t569;
t561 = t442 * t593;
t443 = t489 * t571 + (-pkin(2) * t598 + t489 * t635) * t526 - t492 * t568;
t560 = t443 * t590;
t444 = -t492 * t571 + (-pkin(2) * t604 - t492 * t635) * t526 - t489 * t568;
t559 = t444 * t590;
t481 = mrSges(3,1) * t634;
t507 = mrSges(2,2) - mrSges(3,3);
t533 = (-(t511 * mrSges(3,2) - mrSges(2,1)) * t521 - t507 * t512) * pkin(1);
t448 = (t481 + t576) * t520 + t533 + t567;
t558 = t448 * t587;
t482 = mrSges(3,1) * t633;
t532 = (-(t514 * mrSges(3,2) - mrSges(2,1)) * t524 - t507 * t515) * pkin(1);
t449 = (t482 + t575) * t523 + t532 + t566;
t557 = t449 * t585;
t483 = mrSges(3,1) * t632;
t531 = (-(t517 * mrSges(3,2) - mrSges(2,1)) * t527 - t507 * t518) * pkin(1);
t450 = (t483 + t574) * t526 + t531 + t565;
t556 = t450 * t583;
t555 = t451 * t597;
t554 = t452 * t597;
t553 = t453 * t594;
t552 = t454 * t594;
t551 = t455 * t591;
t550 = t456 * t591;
t549 = t529 * t618;
t548 = t529 * t617;
t547 = t529 * t616;
t463 = t520 * t576 + t567;
t546 = t463 * t587;
t464 = t523 * t575 + t566;
t545 = t464 * t585;
t465 = t526 * t574 + t565;
t544 = t465 * t583;
t472 = Ifges(3,5) * t511 + Ifges(3,6) * t520;
t543 = t472 * t587;
t473 = Ifges(3,5) * t514 + Ifges(3,6) * t523;
t542 = t473 * t585;
t474 = Ifges(3,5) * t517 + Ifges(3,6) * t526;
t541 = t474 * t583;
t540 = t487 * t588;
t539 = t488 * t586;
t538 = t489 * t584;
t537 = t490 * t588;
t536 = t491 * t586;
t535 = t492 * t584;
t534 = Ifges(1,3) + ((m(2) + m(3)) * pkin(1) ^ 2) + t631;
t447 = t480 + 0.2e1 * (t483 + t628) * t526 + 0.2e1 * t531 + t534;
t446 = t479 + 0.2e1 * (t482 + t629) * t523 + 0.2e1 * t532 + t534;
t445 = t478 + 0.2e1 * (t481 + t630) * t520 + 0.2e1 * t533 + t534;
t438 = (t459 * t610 + t474 * t547) * t530;
t437 = (t458 * t611 + t473 * t548) * t530;
t436 = (t457 * t612 + t472 * t549) * t530;
t435 = (t450 * t610 + t465 * t547) * t530;
t434 = (t449 * t611 + t464 * t548) * t530;
t433 = (t448 * t612 + t463 * t549) * t530;
t432 = (t447 * t610 + t450 * t547) * t530;
t431 = (t446 * t611 + t449 * t548) * t530;
t430 = (t445 * t612 + t448 * t549) * t530;
t429 = Ifges(3,3) * t538 + (t444 * t541 + t456 * t619) * t589;
t428 = Ifges(3,3) * t535 + (t443 * t541 + t455 * t619) * t589;
t427 = Ifges(3,3) * t539 + (t442 * t542 + t454 * t620) * t592;
t426 = Ifges(3,3) * t536 + (t441 * t542 + t453 * t620) * t592;
t425 = Ifges(3,3) * t540 + (t440 * t543 + t452 * t621) * t595;
t424 = Ifges(3,3) * t537 + (t439 * t543 + t451 * t621) * t595;
t423 = t474 * t538 + (t444 * t544 + t450 * t622) * t589;
t422 = t474 * t535 + (t443 * t544 + t450 * t623) * t589;
t421 = t473 * t539 + (t442 * t545 + t449 * t624) * t592;
t420 = t473 * t536 + (t441 * t545 + t449 * t625) * t592;
t419 = t472 * t540 + (t440 * t546 + t448 * t626) * t595;
t418 = t472 * t537 + (t439 * t546 + t448 * t627) * t595;
t417 = t459 * t538 + (t444 * t556 + t447 * t622) * t589;
t416 = t459 * t535 + (t443 * t556 + t447 * t623) * t589;
t415 = t458 * t539 + (t442 * t557 + t446 * t624) * t592;
t414 = t458 * t536 + (t441 * t557 + t446 * t625) * t592;
t413 = t457 * t540 + (t440 * t558 + t445 * t626) * t595;
t412 = t457 * t537 + (t439 * t558 + t445 * t627) * t595;
t1 = [m(4) + (t413 * t554 + t415 * t552 + t417 * t550) * t530 + (t425 * t609 + t427 * t607 + t429 * t605 + (t419 * t563 + t421 * t561 + t423 * t559) * t530) * t529, (t413 * t555 + t415 * t553 + t417 * t551) * t530 + (t425 * t603 + t427 * t601 + t429 * t599 + (t419 * t564 + t421 * t562 + t423 * t560) * t530) * t529, (t413 * t612 + t415 * t611 + t417 * t610 + (t419 * t618 + t421 * t617 + t423 * t616) * t529) * t530; (t412 * t554 + t414 * t552 + t416 * t550) * t530 + (t424 * t609 + t426 * t607 + t428 * t605 + (t418 * t563 + t420 * t561 + t422 * t559) * t530) * t529, m(4) + (t412 * t555 + t414 * t553 + t416 * t551) * t530 + (t424 * t603 + t426 * t601 + t428 * t599 + (t418 * t564 + t420 * t562 + t422 * t560) * t530) * t529, (t412 * t612 + t414 * t611 + t416 * t610 + (t418 * t618 + t420 * t617 + t422 * t616) * t529) * t530; (t430 * t554 + t431 * t552 + t432 * t550) * t530 + (t436 * t609 + t437 * t607 + t438 * t605 + (t433 * t563 + t434 * t561 + t435 * t559) * t530) * t529, (t430 * t555 + t431 * t553 + t432 * t551) * t530 + (t436 * t603 + t437 * t601 + t438 * t599 + (t433 * t564 + t434 * t562 + t435 * t560) * t530) * t529, m(4) + (t430 * t612 + t431 * t611 + t432 * t610 + (t433 * t618 + t434 * t617 + t435 * t616) * t529) * t530;];
MX  = t1;
