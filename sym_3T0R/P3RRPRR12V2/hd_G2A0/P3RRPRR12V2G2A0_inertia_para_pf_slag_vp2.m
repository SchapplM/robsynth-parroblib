% Calculate inertia matrix for parallel robot
% P3RRPRR12V2G2A0
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
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V2G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:20:56
% EndTime: 2020-08-06 19:20:58
% DurationCPUTime: 1.91s
% Computational Cost: add. (5925->302), mult. (6831->533), div. (540->6), fcn. (4050->18), ass. (0->190)
t586 = (pkin(1) ^ 2 + pkin(5) ^ 2);
t551 = (pkin(2) + pkin(3));
t635 = -2 * t551;
t537 = sin(qJ(2,3));
t634 = pkin(1) * t537;
t539 = sin(qJ(2,2));
t633 = pkin(1) * t539;
t541 = sin(qJ(2,1));
t632 = pkin(1) * t541;
t631 = mrSges(3,3) - mrSges(2,2);
t630 = (-Ifges(2,1) - Ifges(3,1));
t629 = Ifges(2,6) - Ifges(3,6);
t628 = (mrSges(3,3) * qJ(3,1));
t627 = (mrSges(3,3) * qJ(3,2));
t626 = (mrSges(3,3) * qJ(3,3));
t538 = sin(qJ(1,3));
t550 = pkin(5) - pkin(6);
t500 = t538 * t550;
t543 = cos(qJ(2,3));
t544 = cos(qJ(1,3));
t507 = t537 * qJ(3,3);
t589 = t551 * t543;
t567 = t507 + pkin(1) + t589;
t625 = (t567 * t544 + t500) * t543;
t540 = sin(qJ(1,2));
t501 = t540 * t550;
t545 = cos(qJ(2,2));
t546 = cos(qJ(1,2));
t508 = t539 * qJ(3,2);
t588 = t551 * t545;
t566 = t508 + pkin(1) + t588;
t624 = (t566 * t546 + t501) * t545;
t542 = sin(qJ(1,1));
t502 = t542 * t550;
t547 = cos(qJ(2,1));
t548 = cos(qJ(1,1));
t509 = t541 * qJ(3,1);
t587 = t551 * t547;
t565 = t509 + pkin(1) + t587;
t623 = (t565 * t548 + t502) * t547;
t570 = pkin(1) * t538 - t550 * t544;
t600 = t538 * qJ(3,3);
t573 = t537 * t600;
t622 = (t570 + 0.2e1 * t573) * t551;
t569 = pkin(1) * t540 - t550 * t546;
t597 = t540 * qJ(3,2);
t572 = t539 * t597;
t621 = (t569 + 0.2e1 * t572) * t551;
t568 = pkin(1) * t542 - t550 * t548;
t594 = t542 * qJ(3,1);
t571 = t541 * t594;
t620 = (t568 + 0.2e1 * t571) * t551;
t534 = legFrame(3,2);
t510 = sin(t534);
t619 = t510 * qJ(3,3);
t618 = t510 * t544;
t535 = legFrame(2,2);
t511 = sin(t535);
t617 = t511 * qJ(3,2);
t616 = t511 * t546;
t536 = legFrame(1,2);
t512 = sin(t536);
t615 = t512 * qJ(3,1);
t614 = t512 * t548;
t513 = cos(t534);
t613 = t513 * qJ(3,3);
t612 = t513 * t544;
t514 = cos(t535);
t611 = t514 * qJ(3,2);
t610 = t514 * t546;
t515 = cos(t536);
t609 = t515 * qJ(3,1);
t608 = t515 * t548;
t516 = m(3) * pkin(5) + mrSges(3,2);
t607 = t516 * t537;
t606 = t516 * t539;
t605 = t516 * t541;
t604 = (qJ(3,3) + t551) * (-qJ(3,3) + t551);
t603 = (qJ(3,2) + t551) * (-qJ(3,2) + t551);
t602 = (qJ(3,1) + t551) * (-qJ(3,1) + t551);
t601 = t537 * t551;
t599 = t538 * t551;
t598 = t539 * t551;
t596 = t540 * t551;
t595 = t541 * t551;
t593 = t542 * t551;
t497 = qJ(3,3) + t634;
t592 = t551 * t497;
t498 = qJ(3,2) + t633;
t591 = t551 * t498;
t499 = qJ(3,1) + t632;
t590 = t551 * t499;
t585 = qJ(3,1) * t635;
t584 = qJ(3,2) * t635;
t583 = qJ(3,3) * t635;
t582 = pkin(2) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t530 = 2 * mrSges(3,1) * pkin(2);
t581 = Ifges(3,2) + Ifges(2,3) + t530;
t580 = m(3) * pkin(2) + mrSges(3,1);
t579 = t544 * t607;
t578 = t546 * t606;
t577 = t548 * t605;
t576 = t538 * t604;
t575 = t540 * t603;
t574 = t542 * t602;
t564 = Ifges(2,2) + Ifges(3,3) + t530 + t630;
t563 = Ifges(1,3) + 2 * (mrSges(3,2) + mrSges(2,3)) * pkin(5) - t630 + t586 * m(2);
t561 = pkin(2) ^ 2;
t559 = 0.1e1 / qJ(3,1);
t558 = qJ(3,1) ^ 2;
t557 = 0.1e1 / qJ(3,2);
t556 = qJ(3,2) ^ 2;
t555 = 0.1e1 / qJ(3,3);
t554 = qJ(3,3) ^ 2;
t533 = t547 ^ 2;
t532 = t545 ^ 2;
t531 = t543 ^ 2;
t529 = 2 * t628;
t528 = 2 * t627;
t527 = 2 * t626;
t506 = mrSges(2,1) + t580;
t505 = qJ(3,1) * m(3) + t631;
t504 = m(3) * qJ(3,2) + t631;
t503 = m(3) * qJ(3,3) + t631;
t495 = 0.2e1 * pkin(1) * t506;
t494 = pkin(2) * mrSges(3,2) + t506 * pkin(5) - Ifges(3,4) - Ifges(2,5);
t493 = 0.1e1 / t565;
t492 = 0.1e1 / t566;
t491 = 0.1e1 / t567;
t490 = pkin(1) * qJ(3,1) - t541 * t602;
t489 = pkin(1) * qJ(3,2) - t539 * t603;
t488 = pkin(1) * qJ(3,3) - t537 * t604;
t487 = m(3) * (t558 + t561) + t529 + t581;
t486 = m(3) * (t556 + t561) + t528 + t581;
t485 = m(3) * (t554 + t561) + t527 + t581;
t484 = t568 + t571;
t482 = t569 + t572;
t480 = t570 + t573;
t478 = t568 * t541 + t594;
t477 = t569 * t539 + t597;
t476 = t570 * t537 + t600;
t472 = (qJ(3,1) * mrSges(3,2) + t505 * pkin(5) + t629) * t547 - t541 * t494;
t471 = (qJ(3,2) * mrSges(3,2) + t504 * pkin(5) + t629) * t545 - t539 * t494;
t470 = (qJ(3,3) * mrSges(3,2) + t503 * pkin(5) + t629) * t543 - t537 * t494;
t469 = t548 * t533 * t602 + ((0.2e1 * t509 + pkin(1)) * t548 + t502) * t587 + qJ(3,1) * (t499 * t548 + t541 * t502);
t468 = t546 * t532 * t603 + ((0.2e1 * t508 + pkin(1)) * t546 + t501) * t588 + qJ(3,2) * (t498 * t546 + t539 * t501);
t467 = t544 * t531 * t604 + ((0.2e1 * t507 + pkin(1)) * t544 + t500) * t589 + qJ(3,3) * (t497 * t544 + t537 * t500);
t466 = (t515 * t593 - t615) * t533 + (t484 * t515 + t512 * t595) * t547 + t512 * t499;
t465 = (-t512 * t593 - t609) * t533 + (-t484 * t512 + t515 * t595) * t547 + t515 * t499;
t464 = (t514 * t596 - t617) * t532 + (t482 * t514 + t511 * t598) * t545 + t511 * t498;
t463 = (-t511 * t596 - t611) * t532 + (-t482 * t511 + t514 * t598) * t545 + t514 * t498;
t462 = (t513 * t599 - t619) * t531 + (t480 * t513 + t510 * t601) * t543 + t510 * t497;
t461 = (-t510 * t599 - t613) * t531 + (-t480 * t510 + t513 * t601) * t543 + t513 * t497;
t460 = ((-t558 + t561) * m(3) - (2 * t628) + t564) * t533 + (0.2e1 * (t580 * qJ(3,1) + t582) * t541 + t495) * t547 + 0.2e1 * t505 * t632 + (t558 + t586) * m(3) + t529 + t563;
t459 = ((-t556 + t561) * m(3) - (2 * t627) + t564) * t532 + (0.2e1 * (t580 * qJ(3,2) + t582) * t539 + t495) * t545 + 0.2e1 * t504 * t633 + (t556 + t586) * m(3) + t528 + t563;
t458 = ((-t554 + t561) * m(3) - (2 * t626) + t564) * t531 + (0.2e1 * (t580 * qJ(3,3) + t582) * t537 + t495) * t543 + 0.2e1 * t503 * t634 + (t554 + t586) * m(3) + t527 + t563;
t457 = (t512 * t585 + t515 * t574) * t533 + (-t512 * t490 + t515 * t620) * t547 + t478 * t609 + t512 * t590;
t456 = (-t512 * t574 + t515 * t585) * t533 + (-t515 * t490 - t512 * t620) * t547 - t478 * t615 + t515 * t590;
t455 = (t511 * t584 + t514 * t575) * t532 + (-t511 * t489 + t514 * t621) * t545 + t477 * t611 + t511 * t591;
t454 = (-t511 * t575 + t514 * t584) * t532 + (-t514 * t489 - t511 * t621) * t545 - t477 * t617 + t514 * t591;
t453 = (t510 * t583 + t513 * t576) * t531 + (-t510 * t488 + t513 * t622) * t543 + t476 * t613 + t510 * t592;
t452 = (-t510 * t576 + t513 * t583) * t531 + (-t513 * t488 - t510 * t622) * t543 - t476 * t619 + t513 * t592;
t451 = (-t542 * t605 + (m(3) * t469 - t580 * t623) * t559) * t493;
t450 = (-t540 * t606 + (m(3) * t468 - t580 * t624) * t557) * t492;
t449 = (-t538 * t607 + (m(3) * t467 - t580 * t625) * t555) * t491;
t448 = (-t472 * t542 + (-t469 * t580 + t487 * t623) * t559) * t493;
t447 = (-t471 * t540 + (-t468 * t580 + t486 * t624) * t557) * t492;
t446 = (-t470 * t538 + (-t467 * t580 + t485 * t625) * t555) * t491;
t445 = (t515 * t577 + (m(3) * t457 - t466 * t580) * t559) * t493;
t444 = (t514 * t578 + (m(3) * t455 - t464 * t580) * t557) * t492;
t443 = (t513 * t579 + (m(3) * t453 - t462 * t580) * t555) * t491;
t442 = (-t512 * t577 + (m(3) * t456 - t465 * t580) * t559) * t493;
t441 = (-t511 * t578 + (m(3) * t454 - t463 * t580) * t557) * t492;
t440 = (-t510 * t579 + (m(3) * t452 - t461 * t580) * t555) * t491;
t439 = (-t460 * t542 + (t469 * t605 + t472 * t623) * t559) * t493;
t438 = (-t459 * t540 + (t468 * t606 + t471 * t624) * t557) * t492;
t437 = (-t458 * t538 + (t467 * t607 + t470 * t625) * t555) * t491;
t436 = (t472 * t608 + (-t457 * t580 + t466 * t487) * t559) * t493;
t435 = (t471 * t610 + (-t455 * t580 + t464 * t486) * t557) * t492;
t434 = (t470 * t612 + (-t453 * t580 + t462 * t485) * t555) * t491;
t433 = (-t472 * t614 + (-t456 * t580 + t465 * t487) * t559) * t493;
t432 = (-t471 * t616 + (-t454 * t580 + t463 * t486) * t557) * t492;
t431 = (-t470 * t618 + (-t452 * t580 + t461 * t485) * t555) * t491;
t430 = (t460 * t608 + (t457 * t605 + t466 * t472) * t559) * t493;
t429 = (t459 * t610 + (t455 * t606 + t464 * t471) * t557) * t492;
t428 = (t458 * t612 + (t453 * t607 + t462 * t470) * t555) * t491;
t427 = (-t460 * t614 + (t456 * t605 + t465 * t472) * t559) * t493;
t426 = (-t459 * t616 + (t454 * t606 + t463 * t471) * t557) * t492;
t425 = (-t458 * t618 + (t452 * t607 + t461 * t470) * t555) * t491;
t1 = [m(4) + (t430 * t608 + (t436 * t466 + t445 * t457) * t559) * t493 + (t429 * t610 + (t435 * t464 + t444 * t455) * t557) * t492 + (t428 * t612 + (t434 * t462 + t443 * t453) * t555) * t491, (-t430 * t614 + (t436 * t465 + t445 * t456) * t559) * t493 + (-t429 * t616 + (t435 * t463 + t444 * t454) * t557) * t492 + (-t428 * t618 + (t434 * t461 + t443 * t452) * t555) * t491, (-t430 * t542 + (t436 * t623 + t445 * t469) * t559) * t493 + (-t429 * t540 + (t435 * t624 + t444 * t468) * t557) * t492 + (-t428 * t538 + (t434 * t625 + t443 * t467) * t555) * t491; (t427 * t608 + (t433 * t466 + t442 * t457) * t559) * t493 + (t426 * t610 + (t432 * t464 + t441 * t455) * t557) * t492 + (t425 * t612 + (t431 * t462 + t440 * t453) * t555) * t491, m(4) + (-t427 * t614 + (t433 * t465 + t442 * t456) * t559) * t493 + (-t426 * t616 + (t432 * t463 + t441 * t454) * t557) * t492 + (-t425 * t618 + (t431 * t461 + t440 * t452) * t555) * t491, (-t427 * t542 + (t433 * t623 + t442 * t469) * t559) * t493 + (-t426 * t540 + (t432 * t624 + t441 * t468) * t557) * t492 + (-t425 * t538 + (t431 * t625 + t440 * t467) * t555) * t491; (t439 * t608 + (t448 * t466 + t451 * t457) * t559) * t493 + (t438 * t610 + (t447 * t464 + t450 * t455) * t557) * t492 + (t437 * t612 + (t446 * t462 + t449 * t453) * t555) * t491, (-t439 * t614 + (t448 * t465 + t451 * t456) * t559) * t493 + (-t438 * t616 + (t447 * t463 + t450 * t454) * t557) * t492 + (-t437 * t618 + (t446 * t461 + t449 * t452) * t555) * t491, m(4) + (-t439 * t542 + (t448 * t623 + t451 * t469) * t559) * t493 + (-t438 * t540 + (t447 * t624 + t450 * t468) * t557) * t492 + (-t437 * t538 + (t446 * t625 + t449 * t467) * t555) * t491;];
MX  = t1;
