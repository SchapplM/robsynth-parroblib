% Calculate inertia matrix for parallel robot
% P3RRPRR12V2G3A0
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V2G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:27:23
% EndTime: 2020-08-06 19:27:25
% DurationCPUTime: 2.05s
% Computational Cost: add. (5925->302), mult. (6831->536), div. (540->6), fcn. (4050->18), ass. (0->190)
t586 = (pkin(1) ^ 2 + pkin(5) ^ 2);
t554 = (pkin(2) + pkin(3));
t641 = -2 * t554;
t540 = sin(qJ(2,3));
t640 = pkin(1) * t540;
t542 = sin(qJ(2,2));
t639 = pkin(1) * t542;
t544 = sin(qJ(2,1));
t638 = pkin(1) * t544;
t637 = mrSges(3,3) - mrSges(2,2);
t636 = (-Ifges(2,1) - Ifges(3,1));
t635 = Ifges(2,6) - Ifges(3,6);
t634 = (mrSges(3,3) * qJ(3,1));
t633 = (mrSges(3,3) * qJ(3,2));
t632 = (mrSges(3,3) * qJ(3,3));
t541 = sin(qJ(1,3));
t546 = cos(qJ(2,3));
t510 = t540 * qJ(3,3);
t592 = t554 * t546;
t570 = t510 + pkin(1) + t592;
t547 = cos(qJ(1,3));
t553 = pkin(5) - pkin(6);
t598 = t553 * t547;
t631 = (t570 * t541 - t598) * t546;
t543 = sin(qJ(1,2));
t548 = cos(qJ(2,2));
t511 = t542 * qJ(3,2);
t591 = t554 * t548;
t569 = t511 + pkin(1) + t591;
t549 = cos(qJ(1,2));
t597 = t553 * t549;
t630 = (t569 * t543 - t597) * t548;
t545 = sin(qJ(1,1));
t550 = cos(qJ(2,1));
t512 = t544 * qJ(3,1);
t590 = t554 * t550;
t568 = t512 + pkin(1) + t590;
t551 = cos(qJ(1,1));
t596 = t553 * t551;
t629 = (t568 * t545 - t596) * t550;
t607 = t540 * t547;
t577 = qJ(3,3) * t607;
t589 = pkin(1) * t547 + t541 * t553;
t628 = (0.2e1 * t577 + t589) * t554;
t605 = t542 * t549;
t578 = qJ(3,2) * t605;
t588 = pkin(1) * t549 + t543 * t553;
t627 = (0.2e1 * t578 + t588) * t554;
t603 = t544 * t551;
t579 = qJ(3,1) * t603;
t587 = pkin(1) * t551 + t545 * t553;
t626 = (0.2e1 * t579 + t587) * t554;
t537 = legFrame(3,2);
t513 = sin(t537);
t625 = t513 * qJ(3,3);
t624 = t513 * t541;
t538 = legFrame(2,2);
t514 = sin(t538);
t623 = t514 * qJ(3,2);
t622 = t514 * t543;
t539 = legFrame(1,2);
t515 = sin(t539);
t621 = t515 * qJ(3,1);
t620 = t515 * t545;
t516 = cos(t537);
t619 = t516 * qJ(3,3);
t618 = t516 * t541;
t517 = cos(t538);
t617 = t517 * qJ(3,2);
t616 = t517 * t543;
t518 = cos(t539);
t615 = t518 * qJ(3,1);
t614 = t518 * t545;
t519 = m(3) * pkin(5) + mrSges(3,2);
t613 = t519 * t540;
t612 = t519 * t542;
t611 = t519 * t544;
t610 = (qJ(3,3) + t554) * (-qJ(3,3) + t554);
t609 = (qJ(3,2) + t554) * (-qJ(3,2) + t554);
t608 = (qJ(3,1) + t554) * (-qJ(3,1) + t554);
t606 = t540 * t554;
t604 = t542 * t554;
t602 = t544 * t554;
t601 = t547 * t554;
t600 = t549 * t554;
t599 = t551 * t554;
t500 = qJ(3,3) + t640;
t595 = t554 * t500;
t501 = qJ(3,2) + t639;
t594 = t554 * t501;
t502 = qJ(3,1) + t638;
t593 = t554 * t502;
t585 = qJ(3,1) * t641;
t584 = qJ(3,2) * t641;
t583 = qJ(3,3) * t641;
t582 = pkin(2) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t533 = 2 * mrSges(3,1) * pkin(2);
t581 = Ifges(3,2) + Ifges(2,3) + t533;
t580 = m(3) * pkin(2) + mrSges(3,1);
t576 = t541 * t613;
t575 = t543 * t612;
t574 = t545 * t611;
t573 = t547 * t610;
t572 = t549 * t609;
t571 = t551 * t608;
t567 = Ifges(2,2) + Ifges(3,3) + t533 + t636;
t566 = Ifges(1,3) + 2 * (mrSges(3,2) + mrSges(2,3)) * pkin(5) - t636 + t586 * m(2);
t564 = pkin(2) ^ 2;
t562 = 0.1e1 / qJ(3,1);
t561 = qJ(3,1) ^ 2;
t560 = 0.1e1 / qJ(3,2);
t559 = qJ(3,2) ^ 2;
t558 = 0.1e1 / qJ(3,3);
t557 = qJ(3,3) ^ 2;
t536 = t550 ^ 2;
t535 = t548 ^ 2;
t534 = t546 ^ 2;
t532 = 2 * t634;
t531 = 2 * t633;
t530 = 2 * t632;
t509 = mrSges(2,1) + t580;
t508 = qJ(3,1) * m(3) + t637;
t507 = m(3) * qJ(3,2) + t637;
t506 = m(3) * qJ(3,3) + t637;
t498 = 0.2e1 * pkin(1) * t509;
t497 = pkin(2) * mrSges(3,2) + t509 * pkin(5) - Ifges(3,4) - Ifges(2,5);
t496 = 0.1e1 / t568;
t495 = 0.1e1 / t569;
t494 = 0.1e1 / t570;
t493 = pkin(1) * qJ(3,1) - t544 * t608;
t492 = pkin(1) * qJ(3,2) - t542 * t609;
t491 = pkin(1) * qJ(3,3) - t540 * t610;
t490 = m(3) * (t561 + t564) + t532 + t581;
t489 = m(3) * (t559 + t564) + t531 + t581;
t488 = m(3) * (t557 + t564) + t530 + t581;
t487 = t579 + t587;
t485 = t578 + t588;
t483 = t577 + t589;
t481 = qJ(3,1) * t551 + t587 * t544;
t480 = qJ(3,2) * t549 + t588 * t542;
t479 = qJ(3,3) * t547 + t589 * t540;
t475 = t550 * (qJ(3,1) * mrSges(3,2) + t508 * pkin(5) + t635) - t544 * t497;
t474 = t548 * (qJ(3,2) * mrSges(3,2) + t507 * pkin(5) + t635) - t542 * t497;
t473 = t546 * (qJ(3,3) * mrSges(3,2) + t506 * pkin(5) + t635) - t540 * t497;
t472 = -t545 * t536 * t608 - ((0.2e1 * t512 + pkin(1)) * t545 - t596) * t590 - qJ(3,1) * (t502 * t545 - t544 * t596);
t471 = -t543 * t535 * t609 - ((0.2e1 * t511 + pkin(1)) * t543 - t597) * t591 - qJ(3,2) * (t501 * t543 - t542 * t597);
t470 = -t541 * t534 * t610 - ((0.2e1 * t510 + pkin(1)) * t541 - t598) * t592 - qJ(3,3) * (t500 * t541 - t540 * t598);
t469 = (t518 * t599 - t621) * t536 + (t487 * t518 + t515 * t602) * t550 + t515 * t502;
t468 = (-t515 * t599 - t615) * t536 + (-t487 * t515 + t518 * t602) * t550 + t518 * t502;
t467 = (t517 * t600 - t623) * t535 + (t485 * t517 + t514 * t604) * t548 + t514 * t501;
t466 = (-t514 * t600 - t617) * t535 + (-t485 * t514 + t517 * t604) * t548 + t517 * t501;
t465 = (t516 * t601 - t625) * t534 + (t483 * t516 + t513 * t606) * t546 + t513 * t500;
t464 = (-t513 * t601 - t619) * t534 + (-t483 * t513 + t516 * t606) * t546 + t516 * t500;
t463 = ((-t561 + t564) * m(3) - (2 * t634) + t567) * t536 + (0.2e1 * (t580 * qJ(3,1) + t582) * t544 + t498) * t550 + 0.2e1 * t508 * t638 + (t561 + t586) * m(3) + t532 + t566;
t462 = ((-t559 + t564) * m(3) - (2 * t633) + t567) * t535 + (0.2e1 * (t580 * qJ(3,2) + t582) * t542 + t498) * t548 + 0.2e1 * t507 * t639 + (t559 + t586) * m(3) + t531 + t566;
t461 = ((-t557 + t564) * m(3) - (2 * t632) + t567) * t534 + (0.2e1 * (t580 * qJ(3,3) + t582) * t540 + t498) * t546 + 0.2e1 * t506 * t640 + (t557 + t586) * m(3) + t530 + t566;
t460 = (t515 * t585 + t518 * t571) * t536 + (-t493 * t515 + t518 * t626) * t550 + t481 * t615 + t515 * t593;
t459 = (-t515 * t571 + t518 * t585) * t536 + (-t493 * t518 - t515 * t626) * t550 - t481 * t621 + t518 * t593;
t458 = (t514 * t584 + t517 * t572) * t535 + (-t492 * t514 + t517 * t627) * t548 + t480 * t617 + t514 * t594;
t457 = (-t514 * t572 + t517 * t584) * t535 + (-t492 * t517 - t514 * t627) * t548 - t480 * t623 + t517 * t594;
t456 = (t513 * t583 + t516 * t573) * t534 + (-t491 * t513 + t516 * t628) * t546 + t479 * t619 + t513 * t595;
t455 = (-t513 * t573 + t516 * t583) * t534 + (-t491 * t516 - t513 * t628) * t546 - t479 * t625 + t516 * t595;
t454 = (-t519 * t603 + (m(3) * t472 + t580 * t629) * t562) * t496;
t453 = (-t519 * t605 + (m(3) * t471 + t580 * t630) * t560) * t495;
t452 = (-t519 * t607 + (m(3) * t470 + t580 * t631) * t558) * t494;
t451 = (-t475 * t551 + (-t472 * t580 - t490 * t629) * t562) * t496;
t450 = (-t474 * t549 + (-t471 * t580 - t489 * t630) * t560) * t495;
t449 = (-t473 * t547 + (-t470 * t580 - t488 * t631) * t558) * t494;
t448 = (-t518 * t574 + (m(3) * t460 - t469 * t580) * t562) * t496;
t447 = (-t517 * t575 + (m(3) * t458 - t467 * t580) * t560) * t495;
t446 = (-t516 * t576 + (m(3) * t456 - t465 * t580) * t558) * t494;
t445 = (t515 * t574 + (m(3) * t459 - t468 * t580) * t562) * t496;
t444 = (t514 * t575 + (m(3) * t457 - t466 * t580) * t560) * t495;
t443 = (t513 * t576 + (m(3) * t455 - t464 * t580) * t558) * t494;
t442 = (-t463 * t551 + (t472 * t611 - t475 * t629) * t562) * t496;
t441 = (-t462 * t549 + (t471 * t612 - t474 * t630) * t560) * t495;
t440 = (-t461 * t547 + (t470 * t613 - t473 * t631) * t558) * t494;
t439 = (-t475 * t614 + (-t460 * t580 + t469 * t490) * t562) * t496;
t438 = (-t474 * t616 + (-t458 * t580 + t467 * t489) * t560) * t495;
t437 = (-t473 * t618 + (-t456 * t580 + t465 * t488) * t558) * t494;
t436 = (t475 * t620 + (-t459 * t580 + t468 * t490) * t562) * t496;
t435 = (t474 * t622 + (-t457 * t580 + t466 * t489) * t560) * t495;
t434 = (t473 * t624 + (-t455 * t580 + t464 * t488) * t558) * t494;
t433 = (-t463 * t614 + (t460 * t611 + t469 * t475) * t562) * t496;
t432 = (-t462 * t616 + (t458 * t612 + t467 * t474) * t560) * t495;
t431 = (-t461 * t618 + (t456 * t613 + t465 * t473) * t558) * t494;
t430 = (t463 * t620 + (t459 * t611 + t468 * t475) * t562) * t496;
t429 = (t462 * t622 + (t457 * t612 + t466 * t474) * t560) * t495;
t428 = (t461 * t624 + (t455 * t613 + t464 * t473) * t558) * t494;
t1 = [m(4) + (-t433 * t614 + (t439 * t469 + t448 * t460) * t562) * t496 + (-t432 * t616 + (t438 * t467 + t447 * t458) * t560) * t495 + (-t431 * t618 + (t437 * t465 + t446 * t456) * t558) * t494, (t433 * t620 + (t439 * t468 + t448 * t459) * t562) * t496 + (t432 * t622 + (t438 * t466 + t447 * t457) * t560) * t495 + (t431 * t624 + (t437 * t464 + t446 * t455) * t558) * t494, (-t433 * t551 + (-t439 * t629 + t448 * t472) * t562) * t496 + (-t432 * t549 + (-t438 * t630 + t447 * t471) * t560) * t495 + (-t431 * t547 + (-t437 * t631 + t446 * t470) * t558) * t494; (-t430 * t614 + (t436 * t469 + t445 * t460) * t562) * t496 + (-t429 * t616 + (t435 * t467 + t444 * t458) * t560) * t495 + (-t428 * t618 + (t434 * t465 + t443 * t456) * t558) * t494, m(4) + (t430 * t620 + (t436 * t468 + t445 * t459) * t562) * t496 + (t429 * t622 + (t435 * t466 + t444 * t457) * t560) * t495 + (t428 * t624 + (t434 * t464 + t443 * t455) * t558) * t494, (-t430 * t551 + (-t436 * t629 + t445 * t472) * t562) * t496 + (-t429 * t549 + (-t435 * t630 + t444 * t471) * t560) * t495 + (-t428 * t547 + (-t434 * t631 + t443 * t470) * t558) * t494; (-t442 * t614 + (t451 * t469 + t454 * t460) * t562) * t496 + (-t441 * t616 + (t450 * t467 + t453 * t458) * t560) * t495 + (-t440 * t618 + (t449 * t465 + t452 * t456) * t558) * t494, (t442 * t620 + (t451 * t468 + t454 * t459) * t562) * t496 + (t441 * t622 + (t450 * t466 + t453 * t457) * t560) * t495 + (t440 * t624 + (t449 * t464 + t452 * t455) * t558) * t494, m(4) + (-t442 * t551 + (-t451 * t629 + t454 * t472) * t562) * t496 + (-t441 * t549 + (-t450 * t630 + t453 * t471) * t560) * t495 + (-t440 * t547 + (-t449 * t631 + t452 * t470) * t558) * t494;];
MX  = t1;
