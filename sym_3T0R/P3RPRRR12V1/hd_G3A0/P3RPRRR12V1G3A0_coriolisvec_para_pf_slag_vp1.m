% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR12V1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:00
% EndTime: 2020-08-06 18:28:02
% DurationCPUTime: 1.94s
% Computational Cost: add. (6708->247), mult. (10149->421), div. (2883->7), fcn. (8442->18), ass. (0->161)
t658 = 2 * m(3);
t563 = cos(qJ(1,3));
t556 = sin(qJ(3,3));
t527 = t556 * pkin(3) + qJ(2,3);
t524 = 0.1e1 / t527;
t553 = legFrame(3,2);
t534 = sin(t553);
t537 = cos(t553);
t557 = sin(qJ(1,3));
t570 = xDP(3);
t571 = xDP(2);
t572 = xDP(1);
t493 = (-t557 * t570 + (-t534 * t571 + t537 * t572) * t563) * t524;
t646 = (-pkin(5) - pkin(6));
t543 = pkin(1) - t646;
t545 = 0.1e1 / t556;
t562 = cos(qJ(3,3));
t548 = t562 ^ 2;
t574 = qJ(2,3) ^ 2;
t583 = pkin(3) ^ 2;
t585 = pkin(1) ^ 2;
t654 = 2 * pkin(1);
t586 = -(pkin(6) ^ 2) + (t646 * t654) - t583 - t585 + ((-2 * pkin(6) - pkin(5)) * pkin(5));
t584 = 0.1e1 / pkin(3);
t505 = (-t534 * t572 - t537 * t571) * t584 * t545;
t638 = pkin(3) * t562;
t609 = t505 * t638;
t625 = t543 * t493;
t629 = qJ(2,3) * t556;
t630 = qJ(2,3) * t505;
t641 = pkin(3) * t505;
t588 = qJ(2,3) * t557 + t543 * t563;
t619 = t562 * qJ(2,3);
t622 = t556 * t562;
t481 = t588 * t537 * t556 + t534 * t619 + (t534 * t622 + (-t548 + 0.1e1) * t537 * t557) * pkin(3);
t484 = (-t588 * t534 + t537 * t638) * t556 + t557 * pkin(3) * (t562 - 0.1e1) * (t562 + 0.1e1) * t534 + t537 * t619;
t508 = t527 * t563 - t543 * t557;
t469 = (t508 * t570 + (t481 * t572 + t484 * t571) * t545) * t524;
t651 = 0.2e1 * t469;
t653 = -0.2e1 * pkin(3);
t460 = (((t562 * t625 - t641) * t556 - t630) * t545 * t641 + t625 * t651 + (t543 * t609 + (t548 * t583 + t629 * t653 - t574 + t586) * t493) * t493) * t524;
t463 = (t651 + 0.2e1 * t609 - t625) * t493 * t524;
t568 = rSges(3,3) + pkin(5);
t635 = (-pkin(1) - t568) * m(3);
t518 = -rSges(3,2) * t635 - Icges(3,6);
t519 = -rSges(3,1) * t635 - Icges(3,5);
t499 = -t556 * t518 + t519 * t562;
t502 = t505 ^ 2;
t512 = -t635 + (pkin(1) - rSges(2,2)) * m(2);
t578 = rSges(3,2) ^ 2;
t580 = rSges(3,1) ^ 2;
t516 = (-t578 + t580) * m(3) - Icges(3,1) + Icges(3,2);
t520 = qJ(2,3) * m(3) + m(2) * (rSges(2,3) + qJ(2,3));
t645 = m(3) * rSges(3,2);
t532 = rSges(3,1) * t645 - Icges(3,4);
t587 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - Icges(2,1) - Icges(3,2) - Icges(1,3);
t596 = -rSges(3,1) * t556 - rSges(3,2) * t562;
t597 = rSges(2,3) ^ 2 + t585 + ((-2 * pkin(1) + rSges(2,2)) * rSges(2,2));
t598 = -t580 - t585 + ((-t654 - t568) * t568);
t604 = t502 * t545 * t562;
t647 = -2 * t532;
t648 = 2 * t516;
t652 = 0.2e1 * rSges(2,3);
t616 = (t516 * t548 + t622 * t647 - (qJ(2,3) * t652 + t574 + t597) * m(2) + (0.2e1 * qJ(2,3) * t596 - t574 + t598) * m(3) + t587) * t463 + t512 * t460 + t499 * t604 + (t518 * t562 + t519 * t556) * t502 + (t520 * t651 + (t622 * t648 + (0.4e1 * t548 - 0.2e1) * t532) * t505 + ((rSges(3,1) * t630 + rSges(3,2) * t469) * t562 + (rSges(3,1) * t469 - rSges(3,2) * t630) * t556) * t658) * t493;
t657 = t563 * t616;
t565 = cos(qJ(1,2));
t558 = sin(qJ(3,2));
t528 = t558 * pkin(3) + qJ(2,2);
t525 = 0.1e1 / t528;
t554 = legFrame(2,2);
t535 = sin(t554);
t538 = cos(t554);
t559 = sin(qJ(1,2));
t494 = (-t559 * t570 + (-t535 * t571 + t538 * t572) * t565) * t525;
t546 = 0.1e1 / t558;
t564 = cos(qJ(3,2));
t549 = t564 ^ 2;
t575 = qJ(2,2) ^ 2;
t506 = (-t535 * t572 - t538 * t571) * t584 * t546;
t637 = pkin(3) * t564;
t608 = t506 * t637;
t624 = t543 * t494;
t631 = qJ(2,2) * t558;
t632 = qJ(2,2) * t506;
t640 = pkin(3) * t506;
t589 = qJ(2,2) * t559 + t543 * t565;
t618 = t564 * qJ(2,2);
t621 = t558 * t564;
t482 = t589 * t538 * t558 + t535 * t618 + (t535 * t621 + (-t549 + 0.1e1) * t538 * t559) * pkin(3);
t485 = (-t589 * t535 + t538 * t637) * t558 + t559 * pkin(3) * (t564 - 0.1e1) * (t564 + 0.1e1) * t535 + t538 * t618;
t509 = t528 * t565 - t543 * t559;
t470 = (t509 * t570 + (t482 * t572 + t485 * t571) * t546) * t525;
t650 = 0.2e1 * t470;
t461 = (((t564 * t624 - t640) * t558 - t632) * t546 * t640 + t624 * t650 + (t543 * t608 + (t549 * t583 + t631 * t653 - t575 + t586) * t494) * t494) * t525;
t464 = (t650 + 0.2e1 * t608 - t624) * t494 * t525;
t500 = -t558 * t518 + t519 * t564;
t503 = t506 ^ 2;
t521 = qJ(2,2) * m(3) + m(2) * (rSges(2,3) + qJ(2,2));
t595 = -rSges(3,1) * t558 - rSges(3,2) * t564;
t603 = t503 * t546 * t564;
t615 = (t516 * t549 + t621 * t647 - (qJ(2,2) * t652 + t575 + t597) * m(2) + (0.2e1 * qJ(2,2) * t595 - t575 + t598) * m(3) + t587) * t464 + t512 * t461 + t500 * t603 + (t518 * t564 + t519 * t558) * t503 + (t521 * t650 + (t621 * t648 + (0.4e1 * t549 - 0.2e1) * t532) * t506 + ((rSges(3,1) * t632 + rSges(3,2) * t470) * t564 + (rSges(3,1) * t470 - rSges(3,2) * t632) * t558) * t658) * t494;
t656 = t565 * t615;
t567 = cos(qJ(1,1));
t560 = sin(qJ(3,1));
t529 = t560 * pkin(3) + qJ(2,1);
t526 = 0.1e1 / t529;
t555 = legFrame(1,2);
t536 = sin(t555);
t539 = cos(t555);
t561 = sin(qJ(1,1));
t495 = (-t561 * t570 + (-t536 * t571 + t539 * t572) * t567) * t526;
t547 = 0.1e1 / t560;
t566 = cos(qJ(3,1));
t550 = t566 ^ 2;
t576 = qJ(2,1) ^ 2;
t507 = (-t536 * t572 - t539 * t571) * t584 * t547;
t636 = pkin(3) * t566;
t607 = t507 * t636;
t623 = t543 * t495;
t633 = qJ(2,1) * t560;
t634 = qJ(2,1) * t507;
t639 = pkin(3) * t507;
t590 = qJ(2,1) * t561 + t543 * t567;
t617 = t566 * qJ(2,1);
t620 = t560 * t566;
t483 = t590 * t539 * t560 + t536 * t617 + (t536 * t620 + (-t550 + 0.1e1) * t539 * t561) * pkin(3);
t486 = (-t590 * t536 + t539 * t636) * t560 + t561 * pkin(3) * (t566 - 0.1e1) * (t566 + 0.1e1) * t536 + t539 * t617;
t510 = t529 * t567 - t543 * t561;
t471 = (t510 * t570 + (t483 * t572 + t486 * t571) * t547) * t526;
t649 = 0.2e1 * t471;
t462 = (((t566 * t623 - t639) * t560 - t634) * t547 * t639 + t623 * t649 + (t543 * t607 + (t550 * t583 + t633 * t653 - t576 + t586) * t495) * t495) * t526;
t465 = (t649 + 0.2e1 * t607 - t623) * t495 * t526;
t501 = -t560 * t518 + t519 * t566;
t504 = t507 ^ 2;
t522 = qJ(2,1) * m(3) + m(2) * (rSges(2,3) + qJ(2,1));
t594 = -rSges(3,1) * t560 - rSges(3,2) * t566;
t602 = t504 * t547 * t566;
t614 = (t516 * t550 + t620 * t647 - (qJ(2,1) * t652 + t576 + t597) * m(2) + (0.2e1 * qJ(2,1) * t594 - t576 + t598) * m(3) + t587) * t465 + t512 * t462 + t501 * t602 + (t518 * t566 + t519 * t560) * t504 + (t522 * t649 + (t620 * t648 + (0.4e1 * t550 - 0.2e1) * t532) * t507 + ((rSges(3,1) * t634 + rSges(3,2) * t471) * t566 + (rSges(3,1) * t471 - rSges(3,2) * t634) * t560) * t658) * t495;
t655 = t567 * t614;
t644 = m(3) * (t562 * rSges(3,1) - t556 * rSges(3,2));
t643 = m(3) * (t564 * rSges(3,1) - t558 * rSges(3,2));
t642 = m(3) * (t566 * rSges(3,1) - t560 * rSges(3,2));
t490 = t493 ^ 2;
t573 = -m(2) - m(3);
t613 = t573 * t460 + t512 * t463 - t604 * t644 - t520 * t490 + t596 * (t490 + t502) * m(3);
t491 = t494 ^ 2;
t612 = t573 * t461 + t512 * t464 - t603 * t643 - t521 * t491 + t595 * (t491 + t503) * m(3);
t492 = t495 ^ 2;
t611 = t573 * t462 + t512 * t465 - t602 * t642 - t522 * t492 + t594 * (t492 + t504) * m(3);
t605 = -t645 / 0.2e1;
t610 = rSges(3,1) * t605 + Icges(3,4) / 0.2e1;
t606 = (m(3) * rSges(3,1)) / 0.2e1;
t601 = t613 * t545;
t600 = t612 * t546;
t599 = t611 * t547;
t511 = (t580 / 0.2e1 - t578 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t523 = -(t578 + t580) * m(3) - Icges(3,3);
t593 = (t460 * t644 - t499 * t463 + 0.2e1 * (t532 * t548 + (qJ(2,3) * t606 + t511 * t556) * t562 + t605 * t629 + t610) * t490 - t523 * t604) * t545;
t592 = (t461 * t643 - t500 * t464 + 0.2e1 * (t532 * t549 + (qJ(2,2) * t606 + t511 * t558) * t564 + t605 * t631 + t610) * t491 - t523 * t603) * t546;
t591 = (t462 * t642 - t501 * t465 + 0.2e1 * (t532 * t550 + (qJ(2,1) * t606 + t511 * t560) * t566 + t605 * t633 + t610) * t492 - t523 * t602) * t547;
t1 = [(t483 * t599 + t539 * t655) * t526 + (t482 * t600 + t538 * t656) * t525 + (t481 * t601 + t537 * t657) * t524 + (t534 * t593 + t535 * t592 + t536 * t591) * t584; (t486 * t599 - t536 * t655) * t526 + (t485 * t600 - t535 * t656) * t525 + (t484 * t601 - t534 * t657) * t524 + (t537 * t593 + t538 * t592 + t539 * t591) * t584; (t611 * t510 - t614 * t561) * t526 + (t612 * t509 - t615 * t559) * t525 + (t613 * t508 - t616 * t557) * t524;];
taucX  = t1;
