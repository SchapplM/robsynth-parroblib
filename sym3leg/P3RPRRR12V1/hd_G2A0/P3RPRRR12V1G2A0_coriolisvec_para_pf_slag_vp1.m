% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR12V1G2A0
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
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:28
% EndTime: 2020-08-06 18:24:30
% DurationCPUTime: 1.82s
% Computational Cost: add. (6708->247), mult. (10020->418), div. (2883->7), fcn. (8313->18), ass. (0->164)
t668 = 2 * m(3);
t567 = sin(qJ(1,3));
t566 = sin(qJ(3,3));
t537 = t566 * pkin(3) + qJ(2,3);
t534 = 0.1e1 / t537;
t563 = legFrame(3,2);
t544 = sin(t563);
t547 = cos(t563);
t573 = cos(qJ(1,3));
t580 = xDP(3);
t581 = xDP(2);
t582 = xDP(1);
t500 = (t573 * t580 + (-t544 * t581 + t547 * t582) * t567) * t534;
t656 = (-pkin(5) - pkin(6));
t553 = pkin(1) - t656;
t555 = 0.1e1 / t566;
t572 = cos(qJ(3,3));
t558 = t572 ^ 2;
t584 = qJ(2,3) ^ 2;
t593 = pkin(3) ^ 2;
t595 = pkin(1) ^ 2;
t664 = 2 * pkin(1);
t596 = -(pkin(6) ^ 2) + (t656 * t664) - t593 - t595 + ((-2 * pkin(6) - pkin(5)) * pkin(5));
t594 = 0.1e1 / pkin(3);
t512 = (-t544 * t582 - t547 * t581) * t594 * t555;
t648 = pkin(3) * t572;
t616 = t512 * t648;
t632 = t553 * t500;
t636 = qJ(2,3) * t566;
t637 = qJ(2,3) * t512;
t651 = pkin(3) * t512;
t518 = qJ(2,3) * t573 - t553 * t567;
t626 = t572 * qJ(2,3);
t647 = pkin(3) * t573;
t488 = (-t518 * t547 + t544 * t648) * t566 + (t572 - 0.1e1) * (t572 + 0.1e1) * t547 * t647 + t544 * t626;
t491 = (t518 * t544 + t547 * t648) * t566 + (-t558 + 0.1e1) * t544 * t647 + t547 * t626;
t515 = t567 * t537 + t553 * t573;
t476 = (t515 * t580 + (t488 * t582 + t491 * t581) * t555) * t534;
t661 = 0.2e1 * t476;
t663 = -0.2e1 * pkin(3);
t467 = (((t572 * t632 - t651) * t566 - t637) * t555 * t651 + t632 * t661 + (t553 * t616 + (t558 * t593 + t636 * t663 - t584 + t596) * t500) * t500) * t534;
t470 = (t661 + 0.2e1 * t616 - t632) * t500 * t534;
t578 = rSges(3,3) + pkin(5);
t642 = (-pkin(1) - t578) * m(3);
t528 = -rSges(3,2) * t642 - Icges(3,6);
t529 = -rSges(3,1) * t642 - Icges(3,5);
t506 = -t528 * t566 + t529 * t572;
t509 = t512 ^ 2;
t522 = -t642 + (pkin(1) - rSges(2,2)) * m(2);
t588 = rSges(3,2) ^ 2;
t590 = rSges(3,1) ^ 2;
t526 = (-t588 + t590) * m(3) - Icges(3,1) + Icges(3,2);
t530 = qJ(2,3) * m(3) + m(2) * (rSges(2,3) + qJ(2,3));
t655 = m(3) * rSges(3,2);
t542 = rSges(3,1) * t655 - Icges(3,4);
t597 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - Icges(2,1) - Icges(3,2) - Icges(1,3);
t603 = -rSges(3,1) * t566 - rSges(3,2) * t572;
t604 = rSges(2,3) ^ 2 + t595 + ((-2 * pkin(1) + rSges(2,2)) * rSges(2,2));
t605 = -t590 - t595 + ((-t664 - t578) * t578);
t611 = t509 * t555 * t572;
t629 = t566 * t572;
t657 = -2 * t542;
t658 = 2 * t526;
t662 = 0.2e1 * rSges(2,3);
t623 = (t526 * t558 + t629 * t657 - (qJ(2,3) * t662 + t584 + t604) * m(2) + (0.2e1 * qJ(2,3) * t603 - t584 + t605) * m(3) + t597) * t470 + t522 * t467 + t506 * t611 + (t528 * t572 + t529 * t566) * t509 + (t530 * t661 + (t629 * t658 + (0.4e1 * t558 - 0.2e1) * t542) * t512 + ((rSges(3,1) * t637 + rSges(3,2) * t476) * t572 + (rSges(3,1) * t476 - rSges(3,2) * t637) * t566) * t668) * t500;
t667 = t567 * t623;
t569 = sin(qJ(1,2));
t568 = sin(qJ(3,2));
t538 = t568 * pkin(3) + qJ(2,2);
t535 = 0.1e1 / t538;
t564 = legFrame(2,2);
t545 = sin(t564);
t548 = cos(t564);
t575 = cos(qJ(1,2));
t501 = (t575 * t580 + (-t545 * t581 + t548 * t582) * t569) * t535;
t556 = 0.1e1 / t568;
t574 = cos(qJ(3,2));
t559 = t574 ^ 2;
t585 = qJ(2,2) ^ 2;
t513 = (-t545 * t582 - t548 * t581) * t594 * t556;
t646 = pkin(3) * t574;
t615 = t513 * t646;
t631 = t553 * t501;
t638 = qJ(2,2) * t568;
t639 = qJ(2,2) * t513;
t650 = pkin(3) * t513;
t519 = qJ(2,2) * t575 - t553 * t569;
t625 = t574 * qJ(2,2);
t645 = pkin(3) * t575;
t489 = (-t519 * t548 + t545 * t646) * t568 + (t574 - 0.1e1) * (t574 + 0.1e1) * t548 * t645 + t545 * t625;
t492 = (t519 * t545 + t548 * t646) * t568 + (-t559 + 0.1e1) * t545 * t645 + t548 * t625;
t516 = t569 * t538 + t553 * t575;
t477 = (t516 * t580 + (t489 * t582 + t492 * t581) * t556) * t535;
t660 = 0.2e1 * t477;
t468 = (((t574 * t631 - t650) * t568 - t639) * t556 * t650 + t631 * t660 + (t553 * t615 + (t559 * t593 + t638 * t663 - t585 + t596) * t501) * t501) * t535;
t471 = (t660 + 0.2e1 * t615 - t631) * t501 * t535;
t507 = -t528 * t568 + t529 * t574;
t510 = t513 ^ 2;
t531 = qJ(2,2) * m(3) + m(2) * (rSges(2,3) + qJ(2,2));
t602 = -rSges(3,1) * t568 - rSges(3,2) * t574;
t610 = t510 * t556 * t574;
t628 = t568 * t574;
t622 = (t526 * t559 + t628 * t657 - (qJ(2,2) * t662 + t585 + t604) * m(2) + (0.2e1 * qJ(2,2) * t602 - t585 + t605) * m(3) + t597) * t471 + t522 * t468 + t507 * t610 + (t528 * t574 + t529 * t568) * t510 + (t531 * t660 + (t628 * t658 + (0.4e1 * t559 - 0.2e1) * t542) * t513 + ((rSges(3,1) * t639 + rSges(3,2) * t477) * t574 + (rSges(3,1) * t477 - rSges(3,2) * t639) * t568) * t668) * t501;
t666 = t569 * t622;
t571 = sin(qJ(1,1));
t570 = sin(qJ(3,1));
t539 = t570 * pkin(3) + qJ(2,1);
t536 = 0.1e1 / t539;
t565 = legFrame(1,2);
t546 = sin(t565);
t549 = cos(t565);
t577 = cos(qJ(1,1));
t502 = (t577 * t580 + (-t546 * t581 + t549 * t582) * t571) * t536;
t557 = 0.1e1 / t570;
t576 = cos(qJ(3,1));
t560 = t576 ^ 2;
t586 = qJ(2,1) ^ 2;
t514 = (-t546 * t582 - t549 * t581) * t594 * t557;
t644 = pkin(3) * t576;
t614 = t514 * t644;
t630 = t553 * t502;
t640 = qJ(2,1) * t570;
t641 = qJ(2,1) * t514;
t649 = pkin(3) * t514;
t520 = qJ(2,1) * t577 - t553 * t571;
t624 = t576 * qJ(2,1);
t643 = pkin(3) * t577;
t490 = (-t520 * t549 + t546 * t644) * t570 + (t576 - 0.1e1) * (t576 + 0.1e1) * t549 * t643 + t546 * t624;
t493 = (t520 * t546 + t549 * t644) * t570 + (-t560 + 0.1e1) * t546 * t643 + t549 * t624;
t517 = t571 * t539 + t553 * t577;
t478 = (t517 * t580 + (t490 * t582 + t493 * t581) * t557) * t536;
t659 = 0.2e1 * t478;
t469 = (((t576 * t630 - t649) * t570 - t641) * t557 * t649 + t630 * t659 + (t553 * t614 + (t560 * t593 + t640 * t663 - t586 + t596) * t502) * t502) * t536;
t472 = (t659 + 0.2e1 * t614 - t630) * t502 * t536;
t508 = -t528 * t570 + t529 * t576;
t511 = t514 ^ 2;
t532 = qJ(2,1) * m(3) + m(2) * (rSges(2,3) + qJ(2,1));
t601 = -rSges(3,1) * t570 - rSges(3,2) * t576;
t609 = t511 * t557 * t576;
t627 = t570 * t576;
t621 = (t526 * t560 + t627 * t657 - (qJ(2,1) * t662 + t586 + t604) * m(2) + (0.2e1 * qJ(2,1) * t601 - t586 + t605) * m(3) + t597) * t472 + t522 * t469 + t508 * t609 + (t528 * t576 + t529 * t570) * t511 + (t532 * t659 + (t627 * t658 + (0.4e1 * t560 - 0.2e1) * t542) * t514 + ((rSges(3,1) * t641 + rSges(3,2) * t478) * t576 + (rSges(3,1) * t478 - rSges(3,2) * t641) * t570) * t668) * t502;
t665 = t571 * t621;
t654 = m(3) * (t572 * rSges(3,1) - t566 * rSges(3,2));
t653 = m(3) * (t574 * rSges(3,1) - t568 * rSges(3,2));
t652 = m(3) * (t576 * rSges(3,1) - t570 * rSges(3,2));
t497 = t500 ^ 2;
t583 = -m(2) - m(3);
t620 = t583 * t467 + t522 * t470 - t611 * t654 - t497 * t530 + t603 * (t497 + t509) * m(3);
t498 = t501 ^ 2;
t619 = t583 * t468 + t522 * t471 - t610 * t653 - t498 * t531 + t602 * (t498 + t510) * m(3);
t499 = t502 ^ 2;
t618 = t583 * t469 + t522 * t472 - t609 * t652 - t499 * t532 + t601 * (t499 + t511) * m(3);
t612 = -t655 / 0.2e1;
t617 = rSges(3,1) * t612 + Icges(3,4) / 0.2e1;
t613 = (m(3) * rSges(3,1)) / 0.2e1;
t608 = t620 * t555;
t607 = t619 * t556;
t606 = t618 * t557;
t521 = (t590 / 0.2e1 - t588 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t533 = -(t588 + t590) * m(3) - Icges(3,3);
t600 = (t467 * t654 - t506 * t470 + 0.2e1 * (t542 * t558 + (qJ(2,3) * t613 + t521 * t566) * t572 + t612 * t636 + t617) * t497 - t533 * t611) * t555;
t599 = (t468 * t653 - t507 * t471 + 0.2e1 * (t542 * t559 + (qJ(2,2) * t613 + t521 * t568) * t574 + t612 * t638 + t617) * t498 - t533 * t610) * t556;
t598 = (t469 * t652 - t508 * t472 + 0.2e1 * (t542 * t560 + (qJ(2,1) * t613 + t521 * t570) * t576 + t612 * t640 + t617) * t499 - t533 * t609) * t557;
t1 = [(t490 * t606 + t549 * t665) * t536 + (t489 * t607 + t548 * t666) * t535 + (t488 * t608 + t547 * t667) * t534 + (t544 * t600 + t545 * t599 + t546 * t598) * t594; (t493 * t606 - t546 * t665) * t536 + (t492 * t607 - t545 * t666) * t535 + (t491 * t608 - t544 * t667) * t534 + (t547 * t600 + t548 * t599 + t549 * t598) * t594; (t618 * t517 + t621 * t577) * t536 + (t619 * t516 + t622 * t575) * t535 + (t620 * t515 + t623 * t573) * t534;];
taucX  = t1;
