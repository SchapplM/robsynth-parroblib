% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRR1G2A0
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
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRR1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:24:50
% EndTime: 2020-03-09 21:24:51
% DurationCPUTime: 1.29s
% Computational Cost: add. (14628->201), mult. (9624->301), div. (1818->7), fcn. (6726->59), ass. (0->157)
t633 = 2 * m(3);
t632 = m(3) / 0.2e1;
t568 = 0.1e1 / pkin(3);
t601 = t568 / 0.2e1;
t554 = legFrame(1,2);
t597 = qJ(1,1) + pkin(7);
t525 = t554 + t597;
t517 = qJ(3,1) + t525;
t526 = -t554 + t597;
t518 = qJ(3,1) + t526;
t490 = sin(t517) + sin(t518);
t493 = -cos(t518) + cos(t517);
t549 = pkin(7) + qJ(3,1);
t533 = sin(t549);
t557 = sin(qJ(3,1));
t584 = pkin(1) * t533 + t557 * pkin(2);
t498 = 0.1e1 / t584;
t562 = xDP(3);
t564 = xDP(1);
t602 = t564 / 0.2e1;
t563 = xDP(2);
t603 = t563 / 0.2e1;
t529 = cos(qJ(1,1) + t549);
t604 = t498 * t529;
t454 = t562 * t604 + (t490 * t602 + t493 * t603) * t498;
t591 = t563 * t601;
t541 = qJ(1,1) + t554;
t542 = qJ(1,1) - t554;
t463 = -t493 * pkin(3) + (cos(t526) - cos(t525)) * pkin(2) + (cos(t542) - cos(t541)) * pkin(1);
t613 = t463 * t498;
t457 = t591 * t613;
t600 = t562 * t568;
t487 = -pkin(3) * t529 - pkin(2) * cos(t597) - cos(qJ(1,1)) * pkin(1);
t610 = t487 * t498;
t478 = t600 * t610;
t590 = t564 * t601;
t460 = t490 * pkin(3) + (sin(t525) + sin(t526)) * pkin(2) + (sin(t541) + sin(t542)) * pkin(1);
t616 = t460 * t498;
t578 = t590 * t616;
t442 = -t578 / 0.2e1 + t457 / 0.2e1 + t478 / 0.2e1 + t454;
t448 = t457 + t478 - t578;
t445 = t448 + t454;
t551 = cos(pkin(7));
t623 = pkin(1) * t551;
t530 = pkin(2) + t623;
t536 = cos(t549);
t570 = pkin(1) ^ 2;
t546 = pkin(2) ^ 2 + t570;
t560 = cos(qJ(3,1));
t567 = pkin(3) ^ 2;
t598 = 0.2e1 * pkin(1);
t599 = 0.2e1 * pkin(2) * pkin(3);
t619 = t445 * t448;
t622 = pkin(2) * t551;
t624 = pkin(1) * sin(pkin(7));
t436 = (t442 * t560 * t599 + t445 * t567 + t454 * t546 + (pkin(3) * t442 * t536 + t454 * t622) * t598) * t568 * t498 * t454 + (t530 * t560 - t557 * t624 + pkin(3)) / (t530 * t557 + t560 * t624) * t619;
t581 = -pkin(1) * t536 - pkin(2) * t560;
t439 = (-pkin(3) * t619 + (-pkin(3) * t445 + t581 * t454) * t454) * t498;
t543 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t628 = t581 * rSges(3,1) + t584 * rSges(3,2);
t466 = -Icges(3,3) + (-t543 + t628) * m(3);
t625 = -(t546 + t543) * m(3) - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - Icges(1,3) - Icges(2,3) - Icges(3,3) - 0.2e1 * (m(2) * rSges(2,1) + m(3) * pkin(2)) * t623 + 0.2e1 * m(2) * rSges(2,2) * t624 + (-rSges(2,1) ^ 2 - rSges(2,2) ^ 2 - t570) * m(2);
t430 = t466 * t436 + (t628 * t633 + t625) * t439;
t631 = t498 * t430 / 0.2e1;
t553 = legFrame(2,2);
t596 = qJ(1,2) + pkin(7);
t523 = t553 + t596;
t515 = qJ(3,2) + t523;
t524 = -t553 + t596;
t516 = qJ(3,2) + t524;
t489 = sin(t515) + sin(t516);
t492 = -cos(t516) + cos(t515);
t548 = pkin(7) + qJ(3,2);
t532 = sin(t548);
t556 = sin(qJ(3,2));
t585 = pkin(1) * t532 + t556 * pkin(2);
t497 = 0.1e1 / t585;
t528 = cos(qJ(1,2) + t548);
t606 = t497 * t528;
t453 = t562 * t606 + (t489 * t602 + t492 * t603) * t497;
t539 = qJ(1,2) + t553;
t540 = qJ(1,2) - t553;
t462 = -t492 * pkin(3) + (cos(t524) - cos(t523)) * pkin(2) + (cos(t540) - cos(t539)) * pkin(1);
t614 = t462 * t497;
t456 = t591 * t614;
t486 = -pkin(3) * t528 - pkin(2) * cos(t596) - cos(qJ(1,2)) * pkin(1);
t611 = t486 * t497;
t477 = t600 * t611;
t459 = t489 * pkin(3) + (sin(t523) + sin(t524)) * pkin(2) + (sin(t539) + sin(t540)) * pkin(1);
t617 = t459 * t497;
t579 = t590 * t617;
t441 = -t579 / 0.2e1 + t456 / 0.2e1 + t477 / 0.2e1 + t453;
t447 = t456 + t477 - t579;
t444 = t447 + t453;
t535 = cos(t548);
t559 = cos(qJ(3,2));
t620 = t444 * t447;
t435 = (t441 * t559 * t599 + t444 * t567 + t453 * t546 + (pkin(3) * t441 * t535 + t453 * t622) * t598) * t568 * t497 * t453 + (t530 * t559 - t556 * t624 + pkin(3)) / (t530 * t556 + t559 * t624) * t620;
t582 = -pkin(1) * t535 - pkin(2) * t559;
t437 = (-pkin(3) * t620 + (-pkin(3) * t444 + t582 * t453) * t453) * t497;
t627 = t582 * rSges(3,1) + t585 * rSges(3,2);
t465 = -Icges(3,3) + (-t543 + t627) * m(3);
t429 = t465 * t435 + (t627 * t633 + t625) * t437;
t630 = t497 * t429 / 0.2e1;
t552 = legFrame(3,2);
t595 = qJ(1,3) + pkin(7);
t521 = t552 + t595;
t513 = qJ(3,3) + t521;
t522 = -t552 + t595;
t514 = qJ(3,3) + t522;
t488 = sin(t513) + sin(t514);
t491 = -cos(t514) + cos(t513);
t547 = pkin(7) + qJ(3,3);
t531 = sin(t547);
t555 = sin(qJ(3,3));
t586 = pkin(1) * t531 + t555 * pkin(2);
t496 = 0.1e1 / t586;
t527 = cos(qJ(1,3) + t547);
t608 = t496 * t527;
t452 = t562 * t608 + (t488 * t602 + t491 * t603) * t496;
t537 = qJ(1,3) + t552;
t538 = qJ(1,3) - t552;
t461 = -t491 * pkin(3) + (cos(t522) - cos(t521)) * pkin(2) + (cos(t538) - cos(t537)) * pkin(1);
t615 = t461 * t496;
t455 = t591 * t615;
t485 = -pkin(3) * t527 - pkin(2) * cos(t595) - cos(qJ(1,3)) * pkin(1);
t612 = t485 * t496;
t476 = t600 * t612;
t458 = t488 * pkin(3) + (sin(t521) + sin(t522)) * pkin(2) + (sin(t537) + sin(t538)) * pkin(1);
t618 = t458 * t496;
t580 = t590 * t618;
t440 = -t580 / 0.2e1 + t455 / 0.2e1 + t476 / 0.2e1 + t452;
t446 = t455 + t476 - t580;
t443 = t446 + t452;
t534 = cos(t547);
t558 = cos(qJ(3,3));
t621 = t443 * t446;
t434 = (t440 * t558 * t599 + t443 * t567 + t452 * t546 + (pkin(3) * t440 * t534 + t452 * t622) * t598) * t568 * t496 * t452 + (t530 * t558 - t555 * t624 + pkin(3)) / (t530 * t555 + t558 * t624) * t621;
t583 = -pkin(1) * t534 - pkin(2) * t558;
t438 = (-pkin(3) * t621 + (-pkin(3) * t443 + t583 * t452) * t452) * t496;
t626 = t583 * rSges(3,1) + t586 * rSges(3,2);
t464 = -Icges(3,3) + (-t543 + t626) * m(3);
t428 = t464 * t434 + (t626 * t633 + t625) * t438;
t629 = t496 * t428 / 0.2e1;
t594 = t452 ^ 2 * (pkin(2) * (t555 * rSges(3,1) + rSges(3,2) * t558) + (rSges(3,1) * t531 + rSges(3,2) * t534) * pkin(1)) * t496;
t593 = t453 ^ 2 * (pkin(2) * (t556 * rSges(3,1) + rSges(3,2) * t559) + (rSges(3,1) * t532 + rSges(3,2) * t535) * pkin(1)) * t497;
t592 = t454 ^ 2 * (pkin(2) * (t557 * rSges(3,1) + rSges(3,2) * t560) + (rSges(3,1) * t533 + rSges(3,2) * t536) * pkin(1)) * t498;
t494 = rSges(3,1) * t624 + t530 * rSges(3,2);
t495 = t530 * rSges(3,1) - rSges(3,2) * t624;
t589 = -0.2e1 * t440 * t446 * (t494 * t558 + t555 * t495) * t496;
t588 = -0.2e1 * t441 * t447 * (t494 * t559 + t556 * t495) * t497;
t587 = -0.2e1 * t442 * t448 * (t494 * t560 + t557 * t495) * t498;
t512 = -t543 * m(3) - Icges(3,3);
t433 = t512 * t436 + t466 * t439;
t432 = t512 * t435 + t465 * t437;
t431 = t512 * t434 + t464 * t438;
t1 = [t488 * t629 + t489 * t630 + t490 * t631 + (-t431 * t618 - t432 * t617 - t433 * t616) * t601 + (t488 * t589 + t489 * t588 + t490 * t587 + (-t458 * t594 - t459 * t593 - t460 * t592) * t568) * t632; t491 * t629 + t492 * t630 + t493 * t631 + (t431 * t615 + t432 * t614 + t433 * t613) * t601 + (t491 * t589 + t492 * t588 + t493 * t587 + (t461 * t594 + t462 * t593 + t463 * t592) * t568) * t632; t428 * t608 + t429 * t606 + t430 * t604 + (t431 * t612 + t432 * t611 + t433 * t610) * t568 + (t527 * t589 + t528 * t588 + t529 * t587 + (t485 * t594 + t486 * t593 + t487 * t592) * t568) * m(3);];
taucX  = t1;
