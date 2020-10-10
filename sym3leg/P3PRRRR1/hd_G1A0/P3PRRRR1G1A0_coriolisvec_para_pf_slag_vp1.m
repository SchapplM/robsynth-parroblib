% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR1G1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR1G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:02
% EndTime: 2020-03-09 20:34:03
% DurationCPUTime: 0.85s
% Computational Cost: add. (1206->145), mult. (3351->307), div. (1806->10), fcn. (3585->24), ass. (0->133)
t553 = sin(qJ(3,1));
t559 = cos(qJ(3,1));
t630 = (rSges(3,1) * t559 - rSges(3,2) * t553) * m(3);
t551 = sin(qJ(3,2));
t557 = cos(qJ(3,2));
t629 = (rSges(3,1) * t557 - rSges(3,2) * t551) * m(3);
t549 = sin(qJ(3,3));
t555 = cos(qJ(3,3));
t628 = (rSges(3,1) * t555 - rSges(3,2) * t549) * m(3);
t546 = legFrame(3,3);
t527 = sin(t546);
t530 = cos(t546);
t556 = cos(qJ(2,3));
t606 = t549 * t556;
t501 = -t527 * t606 - t530 * t555;
t505 = -t555 * t527 + t530 * t606;
t564 = xDP(2);
t565 = xDP(1);
t571 = 0.1e1 / pkin(2);
t550 = sin(qJ(2,3));
t534 = 0.1e1 / t550;
t539 = 0.1e1 / t555 ^ 2;
t610 = t534 * t539;
t486 = (t501 * t565 + t505 * t564) * t571 * t610;
t538 = 0.1e1 / t555;
t498 = (t527 * t565 - t530 * t564) * t571 * t538;
t601 = t555 * t556;
t607 = t549 * t555;
t462 = ((-t549 * t550 * t498 + t486 * t601) * t538 * t486 + (-t550 * t486 * t607 + t556 * t498) * t539 * t498) * t534;
t483 = t486 ^ 2;
t495 = t498 ^ 2;
t616 = t495 * t538;
t477 = (-t483 * t555 - t616) * t534 * pkin(2);
t620 = rSges(3,3) * m(3);
t524 = rSges(3,2) * t620 - Icges(3,6);
t525 = rSges(3,1) * t620 - Icges(3,5);
t513 = t524 * t555 + t525 * t549;
t569 = rSges(3,2) ^ 2;
t570 = rSges(3,1) ^ 2;
t595 = t569 + t570;
t522 = -t595 * m(3) - Icges(3,3);
t526 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t537 = t555 ^ 2;
t520 = (-t569 + t570) * m(3) - Icges(3,1) + Icges(3,2);
t590 = t520 * t607;
t593 = t549 * t616;
t516 = t549 * rSges(3,1) + t555 * rSges(3,2);
t613 = t516 * t550;
t621 = 0.2e1 * t526;
t627 = t538 * (m(3) * t477 * t613 + t513 * t462 + t483 * (t537 * t621 - t526 + t590) - t522 * t593);
t547 = legFrame(2,3);
t528 = sin(t547);
t531 = cos(t547);
t558 = cos(qJ(2,2));
t604 = t551 * t558;
t502 = -t528 * t604 - t531 * t557;
t507 = -t557 * t528 + t531 * t604;
t552 = sin(qJ(2,2));
t535 = 0.1e1 / t552;
t542 = 0.1e1 / t557 ^ 2;
t609 = t535 * t542;
t487 = (t502 * t565 + t507 * t564) * t571 * t609;
t541 = 0.1e1 / t557;
t499 = (t528 * t565 - t531 * t564) * t571 * t541;
t600 = t557 * t558;
t605 = t551 * t557;
t463 = ((-t551 * t552 * t499 + t487 * t600) * t541 * t487 + (-t552 * t487 * t605 + t558 * t499) * t542 * t499) * t535;
t484 = t487 ^ 2;
t496 = t499 ^ 2;
t615 = t496 * t541;
t478 = (-t484 * t557 - t615) * t535 * pkin(2);
t514 = t524 * t557 + t525 * t551;
t540 = t557 ^ 2;
t589 = t520 * t605;
t592 = t551 * t615;
t517 = t551 * rSges(3,1) + t557 * rSges(3,2);
t612 = t517 * t552;
t626 = t541 * (m(3) * t478 * t612 + t514 * t463 + t484 * (t540 * t621 - t526 + t589) - t522 * t592);
t548 = legFrame(1,3);
t529 = sin(t548);
t532 = cos(t548);
t560 = cos(qJ(2,1));
t602 = t553 * t560;
t503 = -t529 * t602 - t532 * t559;
t509 = -t559 * t529 + t532 * t602;
t554 = sin(qJ(2,1));
t536 = 0.1e1 / t554;
t545 = 0.1e1 / t559 ^ 2;
t608 = t536 * t545;
t488 = (t503 * t565 + t509 * t564) * t571 * t608;
t544 = 0.1e1 / t559;
t500 = (t529 * t565 - t532 * t564) * t571 * t544;
t599 = t559 * t560;
t603 = t553 * t559;
t464 = ((-t553 * t554 * t500 + t488 * t599) * t544 * t488 + (-t554 * t488 * t603 + t560 * t500) * t545 * t500) * t536;
t485 = t488 ^ 2;
t497 = t500 ^ 2;
t614 = t497 * t544;
t479 = (-t485 * t559 - t614) * t536 * pkin(2);
t515 = t524 * t559 + t525 * t553;
t543 = t559 ^ 2;
t588 = t520 * t603;
t591 = t553 * t614;
t518 = t553 * rSges(3,1) + t559 * rSges(3,2);
t611 = t518 * t554;
t625 = t544 * (m(3) * t479 * t611 + t515 * t464 + t485 * (t543 * t621 - t526 + t588) - t522 * t591);
t523 = m(2) * rSges(2,2) - t620;
t563 = m(2) * rSges(2,1);
t492 = -(t563 + t628) * t556 + t550 * t523;
t533 = -m(1) - m(2) - m(3);
t594 = 0.2e1 * m(3);
t624 = -t556 * (t516 * t498 * t594 + t486 * t523) * t486 + (-t483 * t563 - (t483 + t495) * t628) * t550 + t492 * t462 + t533 * t477;
t493 = -(t563 + t629) * t558 + t552 * t523;
t623 = -t558 * (t517 * t499 * t594 + t487 * t523) * t487 + (-t484 * t563 - (t484 + t496) * t629) * t552 + t493 * t463 + t533 * t478;
t494 = -(t563 + t630) * t560 + t554 * t523;
t622 = -t560 * (t518 * t500 * t594 + t488 * t523) * t488 + (-t485 * t563 - (t485 + t497) * t630) * t554 + t494 * t464 + t533 * t479;
t619 = -t520 / 0.2e1;
t618 = -t524 / 0.4e1;
t617 = t525 / 0.4e1;
t578 = t593 * t613;
t584 = (-m(3) * t578 + t624) * t538 * t534;
t577 = t592 * t612;
t583 = (-m(3) * t577 + t623) * t541 * t535;
t576 = t591 * t611;
t582 = (-m(3) * t576 + t622) * t544 * t536;
t575 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - (0.2e1 * rSges(3,3) ^ 2 + t595) * m(3) / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,3);
t566 = 0.2e1 * qJ(3,3);
t574 = (-0.4e1 * ((t549 * t618 + t555 * t617) * t498 + (t590 / 0.2e1 + (t537 - 0.1e1 / 0.2e1) * t526) * t486) * t498 + t492 * t477 + (cos(t566) * t619 + t526 * sin(t566) + t575) * t462 - t513 * t593) * t610;
t567 = 0.2e1 * qJ(3,2);
t573 = (-0.4e1 * ((t551 * t618 + t557 * t617) * t499 + (t589 / 0.2e1 + (t540 - 0.1e1 / 0.2e1) * t526) * t487) * t499 + t493 * t478 + (cos(t567) * t619 + t526 * sin(t567) + t575) * t463 - t514 * t592) * t609;
t568 = 0.2e1 * qJ(3,1);
t572 = (-0.4e1 * ((t553 * t618 + t559 * t617) * t500 + (t588 / 0.2e1 + (t543 - 0.1e1 / 0.2e1) * t526) * t488) * t500 + t494 * t479 + (cos(t568) * t619 + t526 * sin(t568) + t575) * t464 - t515 * t591) * t608;
t1 = [(t529 * t553 + t532 * t599) * t582 + (t528 * t551 + t531 * t600) * t583 + (t527 * t549 + t530 * t601) * t584 + (t501 * t574 + t502 * t573 + t503 * t572 + t527 * t627 + t528 * t626 + t529 * t625) * t571; (t529 * t599 - t532 * t553) * t582 + (t528 * t600 - t531 * t551) * t583 + (t527 * t601 - t530 * t549) * t584 + (t505 * t574 + t507 * t573 + t509 * t572 - t530 * t627 - t531 * t626 - t532 * t625) * t571; (-t576 - t577 - t578) * m(3) + t622 + t623 + t624;];
taucX  = t1;
