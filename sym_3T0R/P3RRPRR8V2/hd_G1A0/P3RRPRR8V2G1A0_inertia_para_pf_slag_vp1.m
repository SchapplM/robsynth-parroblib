% Calculate inertia matrix for parallel robot
% P3RRPRR8V2G1A0
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:10:57
% EndTime: 2022-11-07 13:10:58
% DurationCPUTime: 1.20s
% Computational Cost: add. (3573->243), mult. (4818->359), div. (249->6), fcn. (1893->44), ass. (0->175)
t638 = rSges(2,3) + pkin(5);
t637 = 2 * pkin(1);
t636 = m(2) / 0.2e1;
t635 = Icges(2,2) / 0.2e1;
t634 = Icges(3,2) / 0.2e1;
t633 = pkin(3) ^ 2 / 0.2e1;
t594 = pkin(2) ^ 2;
t632 = t594 / 0.2e1;
t631 = t637 / 0.2e1;
t630 = pkin(2) * pkin(3);
t629 = -2 * pkin(1);
t627 = m(2) * rSges(2,1);
t626 = m(2) * rSges(2,2);
t625 = m(3) * rSges(3,1);
t555 = qJ(2,3) + pkin(7);
t527 = sin(t555);
t569 = sin(qJ(2,3));
t502 = t569 * pkin(2) + pkin(3) * t527;
t517 = 0.2e1 * t555;
t513 = sin(t517);
t586 = 0.2e1 * qJ(2,3);
t554 = t586 + pkin(7);
t526 = sin(t554);
t561 = sin(t586);
t624 = t502 * t631 + t513 * t633 + t526 * t630 + t561 * t632;
t557 = qJ(2,2) + pkin(7);
t529 = sin(t557);
t571 = sin(qJ(2,2));
t503 = t571 * pkin(2) + pkin(3) * t529;
t518 = 0.2e1 * t557;
t514 = sin(t518);
t587 = 0.2e1 * qJ(2,2);
t556 = t587 + pkin(7);
t528 = sin(t556);
t562 = sin(t587);
t623 = t503 * t631 + t514 * t633 + t528 * t630 + t562 * t632;
t559 = qJ(2,1) + pkin(7);
t531 = sin(t559);
t573 = sin(qJ(2,1));
t504 = t573 * pkin(2) + pkin(3) * t531;
t519 = 0.2e1 * t559;
t515 = sin(t519);
t588 = 0.2e1 * qJ(2,1);
t558 = t588 + pkin(7);
t530 = sin(t558);
t563 = sin(t588);
t622 = t504 * t631 + t515 * t633 + t530 * t630 + t563 * t632;
t590 = rSges(2,2) ^ 2;
t592 = rSges(2,1) ^ 2;
t621 = m(3) * t632 + (-t590 + t592) * t636 + t635 - Icges(2,1) / 0.2e1;
t589 = rSges(3,2) ^ 2;
t591 = rSges(3,1) ^ 2;
t620 = (-t589 + t591) * m(3) / 0.2e1 - Icges(3,1) / 0.2e1 + t634;
t575 = cos(qJ(2,3));
t544 = t575 * pkin(2);
t521 = t544 + pkin(1);
t532 = cos(t555);
t482 = -t532 * rSges(3,1) + t527 * rSges(3,2) - t521;
t619 = m(3) * t482;
t577 = cos(qJ(2,2));
t545 = t577 * pkin(2);
t522 = t545 + pkin(1);
t533 = cos(t557);
t483 = -t533 * rSges(3,1) + t529 * rSges(3,2) - t522;
t618 = m(3) * t483;
t579 = cos(qJ(2,1));
t546 = t579 * pkin(2);
t523 = t546 + pkin(1);
t534 = cos(t559);
t484 = -t534 * rSges(3,1) + t531 * rSges(3,2) - t523;
t617 = m(3) * t484;
t606 = pkin(5) + qJ(3,3);
t550 = -pkin(6) - t606;
t541 = 0.1e1 / t550;
t616 = m(3) * t541;
t607 = pkin(5) + qJ(3,2);
t551 = -pkin(6) - t607;
t542 = 0.1e1 / t551;
t615 = m(3) * t542;
t605 = qJ(3,1) + pkin(5);
t552 = -pkin(6) - t605;
t543 = 0.1e1 / t552;
t614 = m(3) * t543;
t547 = rSges(3,3) + t606;
t613 = m(3) * t547;
t548 = rSges(3,3) + t607;
t612 = m(3) * t548;
t549 = rSges(3,3) + t605;
t611 = m(3) * t549;
t610 = pkin(3) * t532;
t609 = pkin(3) * t533;
t608 = pkin(3) * t534;
t498 = 0.1e1 / (t544 + t610);
t604 = t498 * t541;
t499 = 0.1e1 / (t545 + t609);
t603 = t499 * t542;
t500 = 0.1e1 / (t546 + t608);
t602 = t500 * t543;
t600 = t590 + t592;
t598 = t638 * t626 - Icges(2,6);
t597 = -t638 * t627 + Icges(2,5);
t585 = 2 * pkin(1) ^ 2;
t596 = t585 / 0.2e1 + t589 / 0.2e1 + t591 / 0.2e1 + t632;
t565 = cos(pkin(7));
t516 = t565 * pkin(2) * t625;
t595 = Icges(1,3) + ((2 * pkin(5) ^ 2) + t585 + ((4 * pkin(5) + 2 * rSges(2,3)) * rSges(2,3)) + t600) * t636 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t516 + t634 + t635 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t580 = cos(qJ(1,1));
t578 = cos(qJ(1,2));
t576 = cos(qJ(1,3));
t574 = sin(qJ(1,1));
t572 = sin(qJ(1,2));
t570 = sin(qJ(1,3));
t568 = legFrame(1,3);
t567 = legFrame(2,3);
t566 = legFrame(3,3);
t564 = sin(pkin(7));
t540 = cos(t568);
t539 = cos(t567);
t538 = cos(t566);
t537 = sin(t568);
t536 = sin(t567);
t535 = sin(t566);
t525 = -rSges(2,1) * t626 + Icges(2,4);
t524 = -rSges(3,2) * t625 + Icges(3,4);
t520 = pkin(2) * m(3) + t627;
t512 = rSges(3,1) * t611 - Icges(3,5);
t511 = rSges(3,1) * t612 - Icges(3,5);
t510 = rSges(3,1) * t613 - Icges(3,5);
t509 = rSges(3,2) * t611 - Icges(3,6);
t508 = rSges(3,2) * t612 - Icges(3,6);
t507 = rSges(3,2) * t613 - Icges(3,6);
t496 = -t537 * t574 + t540 * t580;
t495 = t537 * t580 + t540 * t574;
t494 = -t536 * t572 + t539 * t578;
t493 = t536 * t578 + t539 * t572;
t492 = -t535 * t570 + t538 * t576;
t491 = t535 * t576 + t538 * t570;
t490 = t523 * t580 - t574 * t552;
t489 = t522 * t578 - t572 * t551;
t488 = t521 * t576 - t570 * t550;
t487 = t574 * t523 + t580 * t552;
t486 = t572 * t522 + t578 * t551;
t485 = t570 * t521 + t576 * t550;
t481 = 0.2e1 * t516 + t600 * m(2) + Icges(2,3) + Icges(3,3) + (-0.2e1 * t564 * rSges(3,2) * pkin(2) + t589 + t591 + t594) * m(3);
t477 = -t537 * t487 + t490 * t540 + t496 * t608;
t476 = t487 * t540 + t490 * t537 + t495 * t608;
t475 = -t536 * t486 + t489 * t539 + t494 * t609;
t474 = t486 * t539 + t489 * t536 + t493 * t609;
t473 = -t535 * t485 + t488 * t538 + t492 * t610;
t472 = t485 * t538 + t488 * t535 + t491 * t610;
t471 = (-pkin(2) * t611 + t509 * t564 - t512 * t565 + t597) * t573 - (t509 * t565 + t512 * t564 + t598) * t579;
t470 = (-pkin(2) * t612 + t508 * t564 - t511 * t565 + t597) * t571 - (t508 * t565 + t511 * t564 + t598) * t577;
t469 = (-pkin(2) * t613 + t507 * t564 - t510 * t565 + t597) * t569 - (t507 * t565 + t510 * t564 + t598) * t575;
t468 = (t504 * t484 + t622) * m(3) * t602;
t467 = (t503 * t483 + t623) * m(3) * t603;
t466 = (t502 * t482 + t624) * m(3) * t604;
t465 = (t484 * t496 + t477) * t614;
t464 = (t484 * t495 + t476) * t614;
t463 = (t483 * t494 + t475) * t615;
t462 = (t483 * t493 + t474) * t615;
t461 = (t482 * t492 + t473) * t616;
t460 = (t482 * t491 + t472) * t616;
t459 = cos(t519) * t620 + cos(t588) * t621 + (t520 * t579 - t573 * t626) * t637 + t524 * t515 + t525 * t563 + (t549 ^ 2 + (pkin(2) * cos(t558) + t534 * t637) * rSges(3,1) + (t531 * t629 + (-t530 - t564) * pkin(2)) * rSges(3,2) + t596) * m(3) + t595;
t458 = cos(t518) * t620 + cos(t587) * t621 + (t520 * t577 - t571 * t626) * t637 + t524 * t514 + t525 * t562 + (t548 ^ 2 + (pkin(2) * cos(t556) + t533 * t637) * rSges(3,1) + (t529 * t629 + (-t528 - t564) * pkin(2)) * rSges(3,2) + t596) * m(3) + t595;
t457 = cos(t517) * t620 + cos(t586) * t621 + (t520 * t575 - t569 * t626) * t637 + t524 * t513 + t525 * t561 + (t547 ^ 2 + (pkin(2) * cos(t554) + t532 * t637) * rSges(3,1) + (t527 * t629 + (-t526 - t564) * pkin(2)) * rSges(3,2) + t596) * m(3) + t595;
t456 = (t459 * t496 + t477 * t617) * t543;
t455 = (t459 * t495 + t476 * t617) * t543;
t454 = (t458 * t494 + t475 * t618) * t542;
t453 = (t458 * t493 + t474 * t618) * t542;
t452 = (t457 * t492 + t473 * t619) * t541;
t451 = (t457 * t491 + t472 * t619) * t541;
t450 = (t471 - (t504 * t459 + t617 * t622) * t543) * t500;
t449 = (t470 - (t503 * t458 + t618 * t623) * t542) * t499;
t448 = (t469 - (t502 * t457 + t619 * t624) * t541) * t498;
t1 = [m(4) - (-t456 * t496 - t465 * t477) * t543 - (-t454 * t494 - t463 * t475) * t542 - (-t452 * t492 - t461 * t473) * t541, -(-t456 * t495 - t465 * t476) * t543 - (-t454 * t493 - t463 * t474) * t542 - (-t452 * t491 - t461 * t472) * t541, -(-t456 * t504 - t465 * t622 + t496 * t471) * t602 - (-t454 * t503 - t463 * t623 + t494 * t470) * t603 - (-t452 * t502 - t461 * t624 + t492 * t469) * t604; -(-t455 * t496 - t464 * t477) * t543 - (-t453 * t494 - t462 * t475) * t542 - (-t451 * t492 - t460 * t473) * t541, m(4) - (-t455 * t495 - t464 * t476) * t543 - (-t453 * t493 - t462 * t474) * t542 - (-t451 * t491 - t460 * t472) * t541, -(-t455 * t504 - t464 * t622 + t495 * t471) * t602 - (-t453 * t503 - t462 * t623 + t493 * t470) * t603 - (-t451 * t502 - t460 * t624 + t491 * t469) * t604; -(t450 * t496 - t468 * t477) * t543 - (t449 * t494 - t467 * t475) * t542 - (t448 * t492 - t466 * t473) * t541, -(t450 * t495 - t468 * t476) * t543 - (t449 * t493 - t467 * t474) * t542 - (t448 * t491 - t466 * t472) * t541, m(4) + (t500 * t481 - (-t468 * t622 + (t500 * t471 + t450) * t504) * t543) * t500 + (t499 * t481 - (-t467 * t623 + (t499 * t470 + t449) * t503) * t542) * t499 + (t498 * t481 - (-t466 * t624 + (t498 * t469 + t448) * t502) * t541) * t498;];
MX  = t1;
