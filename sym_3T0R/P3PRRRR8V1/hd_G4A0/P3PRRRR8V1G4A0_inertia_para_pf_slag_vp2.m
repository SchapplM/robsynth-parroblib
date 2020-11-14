% Calculate inertia matrix for parallel robot
% P3PRRRR8V1G4A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:27
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V1G4A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:25:12
% EndTime: 2020-08-06 17:25:14
% DurationCPUTime: 1.86s
% Computational Cost: add. (5727->268), mult. (14643->568), div. (648->7), fcn. (16893->34), ass. (0->245)
t665 = 2 * Ifges(3,4);
t568 = sin(qJ(2,3));
t574 = cos(qJ(2,3));
t573 = cos(qJ(3,3));
t658 = pkin(2) * t573;
t512 = -t574 * pkin(5) + t568 * t658;
t553 = sin(pkin(3));
t555 = cos(pkin(3));
t567 = sin(qJ(3,3));
t661 = pkin(2) * t567;
t491 = t512 * t553 + t555 * t661;
t664 = 0.1e1 / t491;
t570 = sin(qJ(2,2));
t576 = cos(qJ(2,2));
t575 = cos(qJ(3,2));
t657 = pkin(2) * t575;
t513 = -t576 * pkin(5) + t570 * t657;
t569 = sin(qJ(3,2));
t660 = pkin(2) * t569;
t492 = t513 * t553 + t555 * t660;
t663 = 0.1e1 / t492;
t572 = sin(qJ(2,1));
t578 = cos(qJ(2,1));
t577 = cos(qJ(3,1));
t656 = pkin(2) * t577;
t514 = -t578 * pkin(5) + t572 * t656;
t571 = sin(qJ(3,1));
t659 = pkin(2) * t571;
t493 = t514 * t553 + t555 * t659;
t662 = 0.1e1 / t493;
t655 = Ifges(3,1) + Ifges(2,3);
t556 = legFrame(3,3);
t527 = sin(t556);
t533 = cos(t556);
t552 = sin(pkin(6));
t554 = cos(pkin(6));
t497 = t527 * t554 + t533 * t552;
t500 = -t527 * t552 + t533 * t554;
t559 = legFrame(3,1);
t530 = sin(t559);
t536 = cos(t559);
t564 = legFrame(3,2);
t542 = sin(t564);
t633 = t536 * t542;
t462 = -t530 * t497 + t500 * t633;
t624 = t555 * t568;
t503 = t552 * t574 + t554 * t624;
t506 = -t552 * t624 + t554 * t574;
t479 = t503 * t533 + t527 * t506;
t585 = t527 * t503 - t506 * t533;
t627 = t553 * t573;
t440 = (t479 * t633 - t585 * t530) * t567 + t462 * t627;
t549 = 0.1e1 / t573;
t654 = t440 * t549;
t557 = legFrame(2,3);
t528 = sin(t557);
t534 = cos(t557);
t498 = t528 * t554 + t534 * t552;
t501 = -t528 * t552 + t534 * t554;
t560 = legFrame(2,1);
t531 = sin(t560);
t537 = cos(t560);
t565 = legFrame(2,2);
t543 = sin(t565);
t632 = t537 * t543;
t464 = -t531 * t498 + t501 * t632;
t623 = t555 * t570;
t504 = t552 * t576 + t554 * t623;
t507 = -t552 * t623 + t554 * t576;
t480 = t504 * t534 + t528 * t507;
t584 = t528 * t504 - t507 * t534;
t626 = t553 * t575;
t441 = (t480 * t632 - t584 * t531) * t569 + t464 * t626;
t550 = 0.1e1 / t575;
t653 = t441 * t550;
t558 = legFrame(1,3);
t529 = sin(t558);
t535 = cos(t558);
t499 = t529 * t554 + t535 * t552;
t502 = -t529 * t552 + t535 * t554;
t561 = legFrame(1,1);
t532 = sin(t561);
t538 = cos(t561);
t566 = legFrame(1,2);
t544 = sin(t566);
t631 = t538 * t544;
t466 = -t532 * t499 + t502 * t631;
t622 = t555 * t572;
t505 = t552 * t578 + t554 * t622;
t508 = -t552 * t622 + t554 * t578;
t481 = t505 * t535 + t529 * t508;
t583 = t529 * t505 - t508 * t535;
t625 = t553 * t577;
t442 = (t481 * t631 - t583 * t532) * t571 + t466 * t625;
t551 = 0.1e1 / t577;
t652 = t442 * t551;
t636 = t530 * t542;
t467 = t497 * t536 + t500 * t636;
t443 = (-t479 * t636 - t585 * t536) * t567 - t467 * t627;
t651 = t443 * t549;
t635 = t531 * t543;
t468 = t498 * t537 + t501 * t635;
t444 = (-t480 * t635 - t584 * t537) * t569 - t468 * t626;
t650 = t444 * t550;
t634 = t532 * t544;
t469 = t499 * t538 + t502 * t634;
t445 = (-t481 * t634 - t583 * t538) * t571 - t469 * t625;
t649 = t445 * t551;
t648 = t664 * t549;
t647 = t663 * t550;
t646 = t662 * t551;
t563 = mrSges(2,2) - mrSges(3,3);
t603 = t573 * mrSges(3,1) - t567 * mrSges(3,2);
t645 = ((mrSges(2,1) + t603) * t574 - t568 * t563) * t553;
t602 = t575 * mrSges(3,1) - t569 * mrSges(3,2);
t644 = ((mrSges(2,1) + t602) * t576 - t570 * t563) * t553;
t601 = t577 * mrSges(3,1) - t571 * mrSges(3,2);
t643 = ((mrSges(2,1) + t601) * t578 - t572 * t563) * t553;
t562 = -Ifges(3,1) + Ifges(3,2);
t509 = (t562 * t573 + t567 * t665) * t573 + t655;
t642 = t509 * t549;
t510 = (t562 * t575 + t569 * t665) * t575 + t655;
t641 = t510 * t550;
t511 = (t562 * t577 + t571 * t665) * t577 + t655;
t640 = t511 * t551;
t518 = Ifges(3,5) * t567 + Ifges(3,6) * t573;
t639 = t518 * t549;
t519 = Ifges(3,5) * t569 + Ifges(3,6) * t575;
t638 = t519 * t550;
t520 = Ifges(3,5) * t571 + Ifges(3,6) * t577;
t637 = t520 * t551;
t545 = cos(t564);
t630 = t545 * t549;
t546 = cos(t565);
t629 = t546 * t550;
t547 = cos(t566);
t628 = t547 * t551;
t621 = t555 * t574;
t620 = t555 * t576;
t619 = t555 * t578;
t461 = t497 * t633 + t530 * t500;
t434 = (-t568 * t461 + t462 * t621) * t658 + pkin(5) * (t574 * t461 + t462 * t624);
t618 = t434 * t648;
t470 = t497 * t636 - t500 * t536;
t435 = -(t467 * t621 - t470 * t568) * t658 - pkin(5) * (t467 * t624 + t574 * t470);
t617 = t435 * t648;
t463 = t498 * t632 + t531 * t501;
t436 = (-t570 * t463 + t464 * t620) * t657 + pkin(5) * (t576 * t463 + t464 * t623);
t616 = t436 * t647;
t471 = t498 * t635 - t501 * t537;
t437 = -(t468 * t620 - t471 * t570) * t657 - pkin(5) * (t468 * t623 + t576 * t471);
t615 = t437 * t647;
t465 = t499 * t631 + t532 * t502;
t438 = (-t572 * t465 + t466 * t619) * t656 + pkin(5) * (t578 * t465 + t466 * t622);
t614 = t438 * t646;
t472 = t499 * t634 - t502 * t538;
t439 = -(t469 * t619 - t472 * t572) * t656 - pkin(5) * (t469 * t622 + t578 * t472);
t613 = t439 * t646;
t458 = t479 * t567 + t500 * t627;
t612 = t458 * t630;
t459 = t480 * t569 + t501 * t626;
t611 = t459 * t629;
t460 = t481 * t571 + t502 * t625;
t610 = t460 * t628;
t579 = 0.1e1 / pkin(2);
t609 = t579 * t648;
t608 = t579 * t647;
t607 = t579 * t646;
t606 = t549 * t645;
t605 = t550 * t644;
t604 = t551 * t643;
t600 = Ifges(3,3) * t609;
t599 = Ifges(3,3) * t608;
t598 = Ifges(3,3) * t607;
t597 = ((-t568 * t497 + t500 * t621) * t658 + pkin(5) * (t574 * t497 + t500 * t624)) * t664 * t630;
t596 = ((-t570 * t498 + t501 * t620) * t657 + pkin(5) * (t576 * t498 + t501 * t623)) * t663 * t629;
t595 = ((-t572 * t499 + t502 * t619) * t656 + pkin(5) * (t578 * t499 + t502 * t622)) * t662 * t628;
t482 = t603 * t555 - (t567 * mrSges(3,1) + t573 * mrSges(3,2)) * t553 * t568;
t594 = t482 * t609;
t483 = t602 * t555 - (t569 * mrSges(3,1) + t575 * mrSges(3,2)) * t553 * t570;
t593 = t483 * t608;
t484 = t601 * t555 - (t571 * mrSges(3,1) + t577 * mrSges(3,2)) * t553 * t572;
t592 = t484 * t607;
t591 = t518 * t609;
t590 = t519 * t608;
t589 = t520 * t607;
t588 = t579 * t597;
t587 = t579 * t596;
t586 = t579 * t595;
t582 = -t512 * t555 + t553 * t661;
t581 = -t513 * t555 + t553 * t660;
t580 = -t514 * t555 + t553 * t659;
t548 = m(1) + m(2) + m(3);
t517 = pkin(5) * t572 + t578 * t656;
t516 = pkin(5) * t570 + t576 * t657;
t515 = pkin(5) * t568 + t574 * t658;
t478 = -t552 * t517 + t580 * t554;
t477 = -t552 * t516 + t581 * t554;
t476 = -t552 * t515 + t582 * t554;
t475 = t517 * t554 + t580 * t552;
t474 = t516 * t554 + t581 * t552;
t473 = t515 * t554 + t582 * t552;
t454 = -t475 * t529 + t478 * t535;
t453 = -t474 * t528 + t477 * t534;
t452 = -t473 * t527 + t476 * t533;
t451 = (t580 * t499 + t502 * t517) * t547 + t544 * t493;
t450 = (t581 * t498 + t501 * t516) * t546 + t543 * t492;
t449 = (t582 * t497 + t500 * t515) * t545 + t542 * t491;
t448 = -t547 * t493 + (t475 * t535 + t478 * t529) * t544;
t447 = -t546 * t492 + (t474 * t534 + t477 * t528) * t543;
t446 = -t545 * t491 + (t473 * t533 + t476 * t527) * t542;
t433 = -t448 * t538 - t532 * t454;
t432 = t448 * t532 - t538 * t454;
t431 = -t447 * t537 - t531 * t453;
t430 = t447 * t531 - t537 * t453;
t429 = -t446 * t536 - t530 * t452;
t428 = t446 * t530 - t536 * t452;
t427 = -Ifges(3,3) * t586 + (t451 * t484 - t520 * t610) * t662;
t426 = -Ifges(3,3) * t587 + (t450 * t483 - t519 * t611) * t663;
t425 = -Ifges(3,3) * t588 + (t449 * t482 - t518 * t612) * t664;
t424 = -t520 * t586 + (t451 * t643 - t511 * t610) * t662;
t423 = -t519 * t587 + (t450 * t644 - t510 * t611) * t663;
t422 = -t518 * t588 + (t449 * t645 - t509 * t612) * t664;
t421 = -t484 * t586 + (-t460 * t547 * t604 + t451 * t548) * t662;
t420 = -t483 * t587 + (-t459 * t546 * t605 + t450 * t548) * t663;
t419 = -t482 * t588 + (-t458 * t545 * t606 + t449 * t548) * t664;
t418 = t438 * t598 + (t433 * t484 + t442 * t637) * t662;
t417 = t439 * t598 + (t432 * t484 + t445 * t637) * t662;
t416 = t436 * t599 + (t431 * t483 + t441 * t638) * t663;
t415 = t437 * t599 + (t430 * t483 + t444 * t638) * t663;
t414 = t434 * t600 + (t429 * t482 + t440 * t639) * t664;
t413 = t435 * t600 + (t428 * t482 + t443 * t639) * t664;
t412 = t438 * t589 + (t433 * t643 + t442 * t640) * t662;
t411 = t439 * t589 + (t432 * t643 + t445 * t640) * t662;
t410 = t436 * t590 + (t431 * t644 + t441 * t641) * t663;
t409 = t437 * t590 + (t430 * t644 + t444 * t641) * t663;
t408 = t434 * t591 + (t429 * t645 + t440 * t642) * t664;
t407 = t435 * t591 + (t428 * t645 + t443 * t642) * t664;
t406 = t438 * t592 + (t433 * t548 + t442 * t604) * t662;
t405 = t439 * t592 + (t432 * t548 + t445 * t604) * t662;
t404 = t436 * t593 + (t431 * t548 + t441 * t605) * t663;
t403 = t437 * t593 + (t430 * t548 + t444 * t605) * t663;
t402 = t434 * t594 + (t429 * t548 + t440 * t606) * t664;
t401 = t435 * t594 + (t428 * t548 + t443 * t606) * t664;
t1 = [m(4) + (t421 * t451 - t424 * t610) * t662 + (t420 * t450 - t423 * t611) * t663 + (t419 * t449 - t422 * t612) * t664 + (-t425 * t597 - t426 * t596 - t427 * t595) * t579, (t421 * t432 + t424 * t649) * t662 + (t420 * t430 + t423 * t650) * t663 + (t419 * t428 + t422 * t651) * t664 + (t425 * t617 + t426 * t615 + t427 * t613) * t579, (t421 * t433 + t424 * t652) * t662 + (t420 * t431 + t423 * t653) * t663 + (t419 * t429 + t422 * t654) * t664 + (t425 * t618 + t426 * t616 + t427 * t614) * t579; (t405 * t451 - t411 * t610) * t662 + (t403 * t450 - t409 * t611) * t663 + (t401 * t449 - t407 * t612) * t664 + (-t413 * t597 - t415 * t596 - t417 * t595) * t579, m(4) + (t405 * t432 + t411 * t649) * t662 + (t403 * t430 + t409 * t650) * t663 + (t401 * t428 + t407 * t651) * t664 + (t413 * t617 + t415 * t615 + t417 * t613) * t579, (t405 * t433 + t411 * t652) * t662 + (t403 * t431 + t409 * t653) * t663 + (t401 * t429 + t407 * t654) * t664 + (t413 * t618 + t415 * t616 + t417 * t614) * t579; (t406 * t451 - t412 * t610) * t662 + (t404 * t450 - t410 * t611) * t663 + (t402 * t449 - t408 * t612) * t664 + (-t414 * t597 - t416 * t596 - t418 * t595) * t579, (t406 * t432 + t412 * t649) * t662 + (t404 * t430 + t410 * t650) * t663 + (t402 * t428 + t408 * t651) * t664 + (t414 * t617 + t416 * t615 + t418 * t613) * t579, m(4) + (t406 * t433 + t412 * t652) * t662 + (t404 * t431 + t410 * t653) * t663 + (t402 * t429 + t408 * t654) * t664 + (t414 * t618 + t416 * t616 + t418 * t614) * t579;];
MX  = t1;
