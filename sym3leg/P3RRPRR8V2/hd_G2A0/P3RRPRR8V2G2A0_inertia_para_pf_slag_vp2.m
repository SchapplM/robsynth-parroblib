% Calculate inertia matrix for parallel robot
% P3RRPRR8V2G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:12:03
% EndTime: 2020-08-06 21:12:04
% DurationCPUTime: 1.34s
% Computational Cost: add. (4653->278), mult. (7920->452), div. (348->9), fcn. (4632->23), ass. (0->197)
t534 = cos(pkin(7));
t528 = t534 ^ 2;
t641 = 0.2e1 * t528;
t640 = 2 * mrSges(3,1);
t639 = 2 * mrSges(3,3);
t548 = cos(qJ(2,3));
t530 = t548 ^ 2;
t638 = 0.2e1 * t530;
t550 = cos(qJ(2,2));
t531 = t550 ^ 2;
t637 = 0.2e1 * t531;
t552 = cos(qJ(2,1));
t532 = t552 ^ 2;
t636 = 0.2e1 * t532;
t635 = 0.2e1 * t534;
t634 = m(3) * pkin(1);
t633 = m(3) * pkin(2);
t542 = sin(qJ(2,3));
t632 = pkin(1) * t542;
t544 = sin(qJ(2,2));
t631 = pkin(1) * t544;
t546 = sin(qJ(2,1));
t630 = pkin(1) * t546;
t543 = sin(qJ(1,3));
t622 = -qJ(3,3) - pkin(5);
t525 = pkin(6) - t622;
t549 = cos(qJ(1,3));
t565 = pkin(1) * t543 - t549 * t525;
t580 = pkin(3) * (t528 - 0.1e1);
t533 = sin(pkin(7));
t596 = t533 * t542;
t629 = pkin(3) * (t543 * t580 + t565 * t596);
t545 = sin(qJ(1,2));
t623 = -qJ(3,2) - pkin(5);
t526 = pkin(6) - t623;
t551 = cos(qJ(1,2));
t564 = pkin(1) * t545 - t551 * t526;
t595 = t533 * t544;
t628 = pkin(3) * (t545 * t580 + t564 * t595);
t547 = sin(qJ(1,1));
t624 = -qJ(3,1) - pkin(5);
t527 = pkin(6) - t624;
t553 = cos(qJ(1,1));
t563 = pkin(1) * t547 - t553 * t527;
t594 = t533 * t546;
t627 = pkin(3) * (t547 * t580 + t563 * t594);
t626 = t533 * pkin(3);
t625 = t534 * pkin(3);
t538 = Ifges(3,2) - Ifges(3,1);
t621 = mrSges(3,2) * t533;
t620 = Ifges(3,4) * t533;
t619 = t533 * mrSges(3,1);
t618 = Ifges(2,5) + (-mrSges(2,1) - t633) * pkin(5);
t502 = pkin(2) + t625;
t583 = pkin(3) * t596;
t477 = 0.1e1 / (t502 * t548 - t583);
t508 = 0.1e1 / t525;
t617 = t477 * t508;
t582 = pkin(3) * t595;
t478 = 0.1e1 / (t502 * t550 - t582);
t509 = 0.1e1 / t526;
t616 = t478 * t509;
t581 = pkin(3) * t594;
t479 = 0.1e1 / (t502 * t552 - t581);
t510 = 0.1e1 / t527;
t615 = t479 * t510;
t569 = pkin(3) * cos(qJ(2,3) + pkin(7)) + t548 * pkin(2);
t483 = 0.1e1 / t569;
t539 = legFrame(3,2);
t511 = sin(t539);
t614 = t483 * t511;
t514 = cos(t539);
t613 = t483 * t514;
t568 = pkin(3) * cos(qJ(2,2) + pkin(7)) + t550 * pkin(2);
t484 = 0.1e1 / t568;
t540 = legFrame(2,2);
t512 = sin(t540);
t612 = t484 * t512;
t515 = cos(t540);
t611 = t484 * t515;
t567 = pkin(3) * cos(qJ(2,1) + pkin(7)) + t552 * pkin(2);
t485 = 0.1e1 / t567;
t541 = legFrame(1,2);
t513 = sin(t541);
t610 = t485 * t513;
t516 = cos(t541);
t609 = t485 * t516;
t608 = t502 * t514;
t607 = t502 * t515;
t606 = t502 * t516;
t605 = t511 * t502;
t604 = t511 * t543;
t603 = t512 * t502;
t602 = t512 * t545;
t601 = t513 * t502;
t600 = t513 * t547;
t599 = t514 * t543;
t598 = t515 * t545;
t597 = t516 * t547;
t499 = t538 * t528;
t557 = pkin(2) ^ 2;
t593 = m(3) * t557 - 0.2e1 * pkin(2) * t621;
t592 = pkin(2) * t640;
t591 = -0.2e1 * (mrSges(2,2) + t619) * pkin(1);
t590 = pkin(2) * t625;
t589 = t514 * t626;
t588 = t515 * t626;
t587 = t516 * t626;
t586 = t511 * t626;
t585 = t512 * t626;
t584 = t513 * t626;
t492 = -t622 * mrSges(3,2) - Ifges(3,6);
t495 = t622 * mrSges(3,1) + Ifges(3,5);
t570 = -pkin(5) * mrSges(2,2) + Ifges(2,6);
t452 = (t492 * t533 + t495 * t534 + (-m(3) * qJ(3,3) - mrSges(3,3)) * pkin(2) + t618) * t542 + t548 * (-t492 * t534 + t495 * t533 + t570);
t579 = t452 * t617;
t493 = -t623 * mrSges(3,2) - Ifges(3,6);
t496 = t623 * mrSges(3,1) + Ifges(3,5);
t453 = (t493 * t533 + t496 * t534 + (-m(3) * qJ(3,2) - mrSges(3,3)) * pkin(2) + t618) * t544 + t550 * (-t493 * t534 + t496 * t533 + t570);
t578 = t453 * t616;
t494 = -t624 * mrSges(3,2) - Ifges(3,6);
t497 = t624 * mrSges(3,1) + Ifges(3,5);
t454 = (t494 * t533 + t497 * t534 + (-m(3) * qJ(3,1) - mrSges(3,3)) * pkin(2) + t618) * t546 + t552 * (-t494 * t534 + t497 * t533 + t570);
t577 = t454 * t615;
t576 = t452 * t614;
t575 = t453 * t612;
t574 = t454 * t610;
t573 = t452 * t613;
t572 = t453 * t611;
t571 = t454 * t609;
t562 = -t621 + t633;
t566 = (t640 * t534 + (2 * mrSges(2,1)) + 0.2e1 * t562) * pkin(1);
t561 = mrSges(3,2) * t534 + t619;
t560 = 0.2e1 * Ifges(3,4) * t641 + 0.2e1 * (-pkin(2) * mrSges(3,2) - t533 * t538) * t534 - 0.2e1 * pkin(2) * t619 + (2 * Ifges(2,4)) - 0.2e1 * Ifges(3,4);
t559 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + 0.2e1 * (mrSges(2,3) + mrSges(3,3)) * pkin(5) + (m(2) + m(3)) * pkin(1) ^ 2 + m(2) * pkin(5) ^ 2 - t499;
t556 = pkin(3) ^ 2;
t558 = t556 * t641 - t556 + t557 + 0.2e1 * t590;
t504 = pkin(1) * t626;
t490 = -t626 + t630;
t489 = -t626 + t631;
t488 = -t626 + t632;
t482 = -mrSges(3,1) * t534 - t562;
t481 = t590 + t557 / 0.2e1 + (t528 - 0.1e1 / 0.2e1) * t556;
t480 = t534 * t592 + Ifges(2,3) + Ifges(3,3) + t593;
t476 = -0.2e1 * t547 * t581 + t563;
t475 = -0.2e1 * t545 * t582 + t564;
t474 = -0.2e1 * t543 * t583 + t565;
t473 = t547 * t527 + (pkin(1) + t567) * t553;
t472 = t545 * t526 + (pkin(1) + t568) * t551;
t471 = t543 * t525 + (pkin(1) + t569) * t549;
t470 = t558 * t546 + t504;
t469 = t558 * t544 + t504;
t468 = t558 * t542 + t504;
t464 = 0.2e1 * t499 + (t592 + 0.4e1 * t620) * t534 + Ifges(2,2) - Ifges(2,1) + t593 - t538;
t463 = t482 * t552 + t561 * t546 - t634;
t462 = t482 * t550 + t561 * t544 - t634;
t461 = t482 * t548 + t561 * t542 - t634;
t460 = (t502 * t597 + t584) * t552 + t546 * (-t547 * t587 + t601);
t459 = (t502 * t598 + t585) * t550 + t544 * (-t545 * t588 + t603);
t458 = (t502 * t599 + t586) * t548 + t542 * (-t543 * t589 + t605);
t457 = (-t502 * t600 + t587) * t552 + (t547 * t584 + t606) * t546;
t456 = (-t502 * t602 + t588) * t550 + (t545 * t585 + t607) * t544;
t455 = (-t502 * t604 + t589) * t548 + (t543 * t586 + t608) * t542;
t451 = (m(3) * t473 + t463 * t553) * t510;
t450 = (m(3) * t472 + t462 * t551) * t509;
t449 = (m(3) * t471 + t461 * t549) * t508;
t448 = t464 * t532 + (t560 * t546 + t566) * t552 + (-mrSges(3,2) * t630 - t620) * t635 + t546 * t591 + t624 ^ 2 * m(3) + qJ(3,1) * t639 + t559;
t447 = t464 * t531 + (t560 * t544 + t566) * t550 + (-mrSges(3,2) * t631 - t620) * t635 + t544 * t591 + t623 ^ 2 * m(3) + qJ(3,2) * t639 + t559;
t446 = t464 * t530 + (t560 * t542 + t566) * t548 + (-mrSges(3,2) * t632 - t620) * t635 + t542 * t591 + t622 ^ 2 * m(3) + qJ(3,3) * t639 + t559;
t445 = (t481 * t597 + t502 * t584) * t636 + (t513 * t470 + t476 * t606) * t552 - t516 * t627 + t490 * t601;
t444 = (-t481 * t600 + t502 * t587) * t636 + (t470 * t516 - t476 * t601) * t552 + t513 * t627 + t490 * t606;
t443 = (t481 * t598 + t502 * t585) * t637 + (t512 * t469 + t475 * t607) * t550 - t515 * t628 + t489 * t603;
t442 = (-t481 * t602 + t502 * t588) * t637 + (t469 * t515 - t475 * t603) * t550 + t512 * t628 + t489 * t607;
t441 = (t481 * t599 + t502 * t586) * t638 + (t511 * t468 + t474 * t608) * t548 - t514 * t629 + t488 * t605;
t440 = (-t481 * t604 + t502 * t589) * t638 + (t468 * t514 - t474 * t605) * t548 + t511 * t629 + t488 * t608;
t439 = t460 * t577 + t480 * t610;
t438 = t459 * t578 + t480 * t612;
t437 = t458 * t579 + t480 * t614;
t436 = t457 * t577 + t480 * t609;
t435 = t456 * t578 + t480 * t611;
t434 = t455 * t579 + t480 * t613;
t433 = (t448 * t553 + t463 * t473) * t510;
t432 = (t447 * t551 + t462 * t472) * t509;
t431 = (t446 * t549 + t461 * t471) * t508;
t430 = (m(3) * t445 + t460 * t463) * t615;
t429 = (m(3) * t443 + t459 * t462) * t616;
t428 = (m(3) * t441 + t458 * t461) * t617;
t427 = (m(3) * t444 + t457 * t463) * t615;
t426 = (m(3) * t442 + t456 * t462) * t616;
t425 = (m(3) * t440 + t455 * t461) * t617;
t424 = t574 + (t445 * t463 + t448 * t460) * t615;
t423 = t575 + (t443 * t462 + t447 * t459) * t616;
t422 = t576 + (t441 * t461 + t446 * t458) * t617;
t421 = t571 + (t444 * t463 + t448 * t457) * t615;
t420 = t572 + (t442 * t462 + t447 * t456) * t616;
t419 = t573 + (t440 * t461 + t446 * t455) * t617;
t1 = [t437 * t614 + t438 * t612 + t439 * t610 + m(4) + (t424 * t460 + t430 * t445) * t615 + (t423 * t459 + t429 * t443) * t616 + (t422 * t458 + t428 * t441) * t617, t437 * t613 + t438 * t611 + t439 * t609 + (t424 * t457 + t430 * t444) * t615 + (t423 * t456 + t429 * t442) * t616 + (t422 * t455 + t428 * t440) * t617, (t424 * t553 + t430 * t473) * t510 + (t423 * t551 + t429 * t472) * t509 + (t422 * t549 + t428 * t471) * t508; t434 * t614 + t435 * t612 + t436 * t610 + (t421 * t460 + t427 * t445) * t615 + (t420 * t459 + t426 * t443) * t616 + (t419 * t458 + t425 * t441) * t617, t434 * t613 + t435 * t611 + t436 * t609 + m(4) + (t421 * t457 + t427 * t444) * t615 + (t420 * t456 + t426 * t442) * t616 + (t419 * t455 + t425 * t440) * t617, (t421 * t553 + t427 * t473) * t510 + (t420 * t551 + t426 * t472) * t509 + (t419 * t549 + t425 * t471) * t508; (t553 * t574 + (t433 * t460 + t445 * t451) * t479) * t510 + (t551 * t575 + (t432 * t459 + t443 * t450) * t478) * t509 + (t549 * t576 + (t431 * t458 + t441 * t449) * t477) * t508, (t553 * t571 + (t433 * t457 + t444 * t451) * t479) * t510 + (t551 * t572 + (t432 * t456 + t442 * t450) * t478) * t509 + (t549 * t573 + (t431 * t455 + t440 * t449) * t477) * t508, m(4) + (t433 * t553 + t451 * t473) * t510 + (t432 * t551 + t450 * t472) * t509 + (t431 * t549 + t449 * t471) * t508;];
MX  = t1;
