% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRR1G3P3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:26:46
% EndTime: 2020-03-09 21:26:48
% DurationCPUTime: 1.42s
% Computational Cost: add. (14628->213), mult. (9624->317), div. (1818->7), fcn. (6726->59), ass. (0->164)
t644 = 2 * m(3);
t643 = m(3) / 0.2e1;
t572 = 0.1e1 / pkin(3);
t615 = t572 / 0.2e1;
t568 = xDP(1);
t598 = t568 * t615;
t558 = legFrame(1,2);
t605 = qJ(1,1) + pkin(7);
t529 = t558 + t605;
t521 = qJ(3,1) + t529;
t530 = -t558 + t605;
t522 = qJ(3,1) + t530;
t497 = cos(t522) + cos(t521);
t545 = qJ(1,1) + t558;
t546 = qJ(1,1) - t558;
t470 = -t497 * pkin(3) + (-cos(t530) - cos(t529)) * pkin(2) + (-cos(t546) - cos(t545)) * pkin(1);
t553 = pkin(7) + qJ(3,1);
t537 = sin(t553);
t561 = sin(qJ(3,1));
t585 = pkin(1) * t537 + t561 * pkin(2);
t502 = 0.1e1 / t585;
t627 = t470 * t502;
t464 = t598 * t627;
t494 = -sin(t521) + sin(t522);
t467 = t494 * pkin(3) + (-sin(t529) + sin(t530)) * pkin(2) + (-sin(t545) + sin(t546)) * pkin(1);
t566 = xDP(3);
t614 = t566 * t572;
t533 = sin(qJ(1,1) + t553);
t491 = pkin(3) * t533 + pkin(2) * sin(t605) + sin(qJ(1,1)) * pkin(1);
t624 = t491 * t502;
t485 = t614 * t624;
t567 = xDP(2);
t599 = t567 * t615;
t591 = -t599 / 0.2e1;
t616 = t568 / 0.2e1;
t617 = t567 / 0.2e1;
t607 = (t494 * t617 + t497 * t616) * t502;
t618 = t533 * t566;
t449 = t464 / 0.2e1 + t485 / 0.2e1 + (t467 * t591 - t618) * t502 + t607;
t592 = t467 * t599;
t610 = t464 + t485;
t452 = (-t592 - t618) * t502 + t607 + t610;
t461 = -t502 * t618 + t607;
t555 = cos(pkin(7));
t634 = pkin(1) * t555;
t534 = pkin(2) + t634;
t540 = cos(t553);
t574 = pkin(1) ^ 2;
t550 = pkin(2) ^ 2 + t574;
t564 = cos(qJ(3,1));
t571 = pkin(3) ^ 2;
t606 = 0.2e1 * pkin(1);
t613 = 0.2e1 * pkin(2) * pkin(3);
t455 = -t502 * t592 + t610;
t630 = t452 * t455;
t633 = pkin(2) * t555;
t635 = pkin(1) * sin(pkin(7));
t443 = (t449 * t564 * t613 + t452 * t571 + t461 * t550 + (pkin(3) * t449 * t540 + t461 * t633) * t606) * t572 * t502 * t461 + (t534 * t564 - t561 * t635 + pkin(3)) / (t534 * t561 + t564 * t635) * t630;
t582 = -pkin(1) * t540 - pkin(2) * t564;
t445 = (-pkin(3) * t630 + (-pkin(3) * t452 + t582 * t461) * t461) * t502;
t547 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t639 = t582 * rSges(3,1) + t585 * rSges(3,2);
t473 = -Icges(3,3) + (-t547 + t639) * m(3);
t636 = -(t550 + t547) * m(3) - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - Icges(1,3) - Icges(2,3) - Icges(3,3) - 0.2e1 * (m(2) * rSges(2,1) + m(3) * pkin(2)) * t634 + 0.2e1 * m(2) * rSges(2,2) * t635 + (-rSges(2,1) ^ 2 - rSges(2,2) ^ 2 - t574) * m(2);
t621 = t502 * (t473 * t443 + (t639 * t644 + t636) * t445);
t642 = t621 / 0.2e1;
t557 = legFrame(2,2);
t604 = qJ(1,2) + pkin(7);
t527 = t557 + t604;
t519 = qJ(3,2) + t527;
t528 = -t557 + t604;
t520 = qJ(3,2) + t528;
t496 = cos(t520) + cos(t519);
t543 = qJ(1,2) + t557;
t544 = qJ(1,2) - t557;
t469 = -t496 * pkin(3) + (-cos(t528) - cos(t527)) * pkin(2) + (-cos(t544) - cos(t543)) * pkin(1);
t552 = pkin(7) + qJ(3,2);
t536 = sin(t552);
t560 = sin(qJ(3,2));
t586 = pkin(1) * t536 + t560 * pkin(2);
t501 = 0.1e1 / t586;
t628 = t469 * t501;
t463 = t598 * t628;
t493 = -sin(t519) + sin(t520);
t466 = t493 * pkin(3) + (-sin(t527) + sin(t528)) * pkin(2) + (-sin(t543) + sin(t544)) * pkin(1);
t532 = sin(qJ(1,2) + t552);
t490 = pkin(3) * t532 + pkin(2) * sin(t604) + sin(qJ(1,2)) * pkin(1);
t625 = t490 * t501;
t484 = t614 * t625;
t608 = (t493 * t617 + t496 * t616) * t501;
t619 = t532 * t566;
t448 = t463 / 0.2e1 + t484 / 0.2e1 + (t466 * t591 - t619) * t501 + t608;
t593 = t466 * t599;
t611 = t463 + t484;
t451 = (-t593 - t619) * t501 + t608 + t611;
t460 = -t501 * t619 + t608;
t539 = cos(t552);
t563 = cos(qJ(3,2));
t454 = -t501 * t593 + t611;
t631 = t451 * t454;
t442 = (t448 * t563 * t613 + t451 * t571 + t460 * t550 + (pkin(3) * t448 * t539 + t460 * t633) * t606) * t572 * t501 * t460 + (t534 * t563 - t560 * t635 + pkin(3)) / (t534 * t560 + t563 * t635) * t631;
t583 = -pkin(1) * t539 - pkin(2) * t563;
t444 = (-pkin(3) * t631 + (-pkin(3) * t451 + t583 * t460) * t460) * t501;
t638 = t583 * rSges(3,1) + t586 * rSges(3,2);
t472 = -Icges(3,3) + (-t547 + t638) * m(3);
t622 = t501 * (t472 * t442 + (t638 * t644 + t636) * t444);
t641 = t622 / 0.2e1;
t556 = legFrame(3,2);
t603 = qJ(1,3) + pkin(7);
t525 = t556 + t603;
t517 = qJ(3,3) + t525;
t526 = -t556 + t603;
t518 = qJ(3,3) + t526;
t495 = cos(t518) + cos(t517);
t541 = qJ(1,3) + t556;
t542 = qJ(1,3) - t556;
t468 = -t495 * pkin(3) + (-cos(t526) - cos(t525)) * pkin(2) + (-cos(t542) - cos(t541)) * pkin(1);
t551 = pkin(7) + qJ(3,3);
t535 = sin(t551);
t559 = sin(qJ(3,3));
t587 = pkin(1) * t535 + t559 * pkin(2);
t500 = 0.1e1 / t587;
t629 = t468 * t500;
t462 = t598 * t629;
t492 = -sin(t517) + sin(t518);
t465 = t492 * pkin(3) + (-sin(t525) + sin(t526)) * pkin(2) + (-sin(t541) + sin(t542)) * pkin(1);
t531 = sin(qJ(1,3) + t551);
t489 = pkin(3) * t531 + pkin(2) * sin(t603) + sin(qJ(1,3)) * pkin(1);
t626 = t489 * t500;
t483 = t614 * t626;
t609 = (t492 * t617 + t495 * t616) * t500;
t620 = t531 * t566;
t447 = t462 / 0.2e1 + t483 / 0.2e1 + (t465 * t591 - t620) * t500 + t609;
t594 = t465 * t599;
t612 = t462 + t483;
t450 = (-t594 - t620) * t500 + t609 + t612;
t459 = -t500 * t620 + t609;
t538 = cos(t551);
t562 = cos(qJ(3,3));
t453 = -t500 * t594 + t612;
t632 = t450 * t453;
t441 = (t447 * t562 * t613 + t450 * t571 + t459 * t550 + (pkin(3) * t447 * t538 + t459 * t633) * t606) * t572 * t500 * t459 + (t534 * t562 - t559 * t635 + pkin(3)) / (t534 * t559 + t562 * t635) * t632;
t584 = -pkin(1) * t538 - pkin(2) * t562;
t446 = (-pkin(3) * t632 + (-pkin(3) * t450 + t584 * t459) * t459) * t500;
t637 = t584 * rSges(3,1) + t587 * rSges(3,2);
t471 = -Icges(3,3) + (-t547 + t637) * m(3);
t623 = t500 * (t471 * t441 + (t637 * t644 + t636) * t446);
t640 = t623 / 0.2e1;
t602 = t459 ^ 2 * (pkin(2) * (t559 * rSges(3,1) + rSges(3,2) * t562) + (rSges(3,1) * t535 + rSges(3,2) * t538) * pkin(1)) * t500;
t601 = t460 ^ 2 * (pkin(2) * (t560 * rSges(3,1) + rSges(3,2) * t563) + (rSges(3,1) * t536 + rSges(3,2) * t539) * pkin(1)) * t501;
t600 = t461 ^ 2 * (pkin(2) * (t561 * rSges(3,1) + rSges(3,2) * t564) + (rSges(3,1) * t537 + rSges(3,2) * t540) * pkin(1)) * t502;
t498 = rSges(3,1) * t635 + t534 * rSges(3,2);
t499 = t534 * rSges(3,1) - rSges(3,2) * t635;
t597 = t447 * t453 * (t498 * t562 + t559 * t499) * t500;
t596 = t448 * t454 * (t498 * t563 + t560 * t499) * t501;
t595 = t449 * t455 * (t498 * t564 + t561 * t499) * t502;
t590 = -0.2e1 * t597;
t589 = -0.2e1 * t596;
t588 = -0.2e1 * t595;
t516 = -t547 * m(3) - Icges(3,3);
t440 = t516 * t443 + t473 * t445;
t439 = t516 * t442 + t472 * t444;
t438 = t516 * t441 + t471 * t446;
t1 = [t495 * t640 + t496 * t641 + t497 * t642 + (t438 * t629 + t439 * t628 + t440 * t627) * t615 + (t495 * t590 + t496 * t589 + t497 * t588 + (t468 * t602 + t469 * t601 + t470 * t600) * t572) * t643; t492 * t640 + t493 * t641 + t494 * t642 + (-t465 * t500 * t438 - t466 * t501 * t439 - t467 * t502 * t440) * t615 + (t492 * t590 + t493 * t589 + t494 * t588 + (-t465 * t602 - t466 * t601 - t467 * t600) * t572) * t643; -t531 * t623 - t532 * t622 - t533 * t621 + (t438 * t626 + t439 * t625 + t440 * t624) * t572 + (0.2e1 * t531 * t597 + 0.2e1 * t532 * t596 + 0.2e1 * t533 * t595 + (t489 * t602 + t490 * t601 + t491 * t600) * t572) * m(3);];
taucX  = t1;
