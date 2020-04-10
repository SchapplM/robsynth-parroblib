% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRR2G2P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR2G2P3A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:40
% EndTime: 2020-03-09 21:21:42
% DurationCPUTime: 1.33s
% Computational Cost: add. (9456->159), mult. (10017->321), div. (2034->9), fcn. (7461->33), ass. (0->191)
t554 = sin(qJ(3,3));
t543 = 0.1e1 / t554 ^ 2;
t560 = cos(qJ(3,3));
t671 = t543 * t560;
t556 = sin(qJ(3,2));
t546 = 0.1e1 / t556 ^ 2;
t562 = cos(qJ(3,2));
t670 = t546 * t562;
t558 = sin(qJ(3,1));
t549 = 0.1e1 / t558 ^ 2;
t564 = cos(qJ(3,1));
t669 = t549 * t564;
t570 = 0.1e1 / pkin(2);
t668 = (pkin(1) * t560 + pkin(2)) * t570;
t667 = (pkin(1) * t562 + pkin(2)) * t570;
t666 = (pkin(1) * t564 + pkin(2)) * t570;
t542 = 0.1e1 / t554;
t665 = t542 * t671;
t545 = 0.1e1 / t556;
t664 = t545 * t670;
t548 = 0.1e1 / t558;
t663 = t548 * t669;
t566 = xDP(3);
t662 = -2 * t566;
t551 = legFrame(3,2);
t621 = qJ(2,3) + qJ(3,3);
t521 = t551 + t621;
t522 = -t551 + t621;
t533 = sin(t621);
t567 = xDP(2);
t568 = xDP(1);
t506 = t533 * t662 + (sin(t521) - sin(t522)) * t568 + (cos(t522) + cos(t521)) * t567;
t661 = -t506 / 0.2e1;
t552 = legFrame(2,2);
t622 = qJ(2,2) + qJ(3,2);
t523 = t552 + t622;
t524 = -t552 + t622;
t534 = sin(t622);
t507 = t534 * t662 + (-sin(t524) + sin(t523)) * t568 + (cos(t524) + cos(t523)) * t567;
t660 = -t507 / 0.2e1;
t553 = legFrame(1,2);
t623 = qJ(2,1) + qJ(3,1);
t525 = t553 + t623;
t526 = -t553 + t623;
t535 = sin(t623);
t508 = t535 * t662 + (sin(t525) - sin(t526)) * t568 + (cos(t526) + cos(t525)) * t567;
t659 = -t508 / 0.2e1;
t572 = 0.1e1 / pkin(1) ^ 2;
t658 = -t572 / 0.4e1;
t657 = pkin(2) * t560;
t656 = pkin(2) * t562;
t655 = pkin(2) * t564;
t654 = pkin(2) * t566;
t633 = t543 * t572;
t536 = sin(t551);
t539 = cos(t551);
t512 = t536 * t568 + t539 * t567;
t530 = pkin(1) + t657;
t555 = sin(qJ(2,3));
t561 = cos(qJ(2,3));
t497 = (-t512 * t530 + t554 * t654) * t561 + t555 * (pkin(2) * t512 * t554 + t530 * t566);
t571 = 0.1e1 / pkin(1);
t634 = t542 * t571;
t614 = t497 * t634;
t494 = t570 * t614;
t596 = t506 * t542 / 0.2e1;
t500 = t571 * t596;
t488 = t500 + t494;
t650 = t488 * t497;
t482 = t633 * t650;
t635 = t542 * t560;
t485 = -pkin(2) * t488 + t635 * t661;
t587 = t633 * t661;
t473 = t485 * t587 + t482;
t569 = pkin(2) ^ 2;
t476 = (t488 * t569 + (0.2e1 * (t500 + t494 / 0.2e1) * t657 + t596) * pkin(1)) * t570 * t587;
t467 = -t482 * t668 + t473 + t476;
t653 = t467 * t542;
t630 = t546 * t572;
t537 = sin(t552);
t540 = cos(t552);
t513 = t537 * t568 + t540 * t567;
t531 = pkin(1) + t656;
t557 = sin(qJ(2,2));
t563 = cos(qJ(2,2));
t498 = (-t513 * t531 + t556 * t654) * t563 + t557 * (pkin(2) * t513 * t556 + t531 * t566);
t631 = t545 * t571;
t613 = t498 * t631;
t495 = t570 * t613;
t595 = t507 * t545 / 0.2e1;
t501 = t571 * t595;
t489 = t501 + t495;
t649 = t489 * t498;
t483 = t630 * t649;
t632 = t545 * t562;
t486 = -pkin(2) * t489 + t632 * t660;
t586 = t630 * t660;
t474 = t486 * t586 + t483;
t477 = (t489 * t569 + (0.2e1 * (t501 + t495 / 0.2e1) * t656 + t595) * pkin(1)) * t570 * t586;
t468 = -t483 * t667 + t474 + t477;
t652 = t468 * t545;
t627 = t549 * t572;
t538 = sin(t553);
t541 = cos(t553);
t514 = t538 * t568 + t541 * t567;
t532 = pkin(1) + t655;
t559 = sin(qJ(2,1));
t565 = cos(qJ(2,1));
t499 = (-t514 * t532 + t558 * t654) * t565 + t559 * (pkin(2) * t514 * t558 + t532 * t566);
t628 = t548 * t571;
t612 = t499 * t628;
t496 = t570 * t612;
t594 = t508 * t548 / 0.2e1;
t502 = t571 * t594;
t490 = t502 + t496;
t648 = t490 * t499;
t484 = t627 * t648;
t629 = t548 * t564;
t487 = -t490 * pkin(2) + t629 * t659;
t585 = t627 * t659;
t475 = t487 * t585 + t484;
t478 = (t490 * t569 + (0.2e1 * (t502 + t496 / 0.2e1) * t655 + t594) * pkin(1)) * t570 * t585;
t469 = -t484 * t666 + t475 + t478;
t651 = t469 * t548;
t491 = t506 * t634 + t494;
t647 = t491 * t497;
t492 = t507 * t631 + t495;
t646 = t492 * t498;
t493 = t508 * t628 + t496;
t645 = t493 * t499;
t470 = t476 + 0.2e1 * t482 + (-t485 * t506 - t650 * t668) * t633;
t626 = t554 * t555;
t515 = t560 * t561 - t626;
t644 = t515 * t470;
t643 = t515 * t542;
t471 = t477 + 0.2e1 * t483 + (-t486 * t507 - t649 * t667) * t630;
t625 = t556 * t557;
t516 = t562 * t563 - t625;
t642 = t516 * t471;
t641 = t516 * t545;
t472 = t478 + 0.2e1 * t484 + (-t487 * t508 - t648 * t666) * t627;
t624 = t558 * t559;
t517 = t564 * t565 - t624;
t640 = t517 * t472;
t639 = t517 * t548;
t638 = t533 * t542;
t637 = t534 * t545;
t636 = t535 * t548;
t509 = -pkin(2) * t626 + t530 * t561;
t620 = t509 * t653;
t510 = -pkin(2) * t625 + t531 * t563;
t619 = t510 * t652;
t511 = -pkin(2) * t624 + t532 * t565;
t618 = t511 * t651;
t617 = t470 * t635;
t616 = t471 * t632;
t615 = t472 * t629;
t611 = t536 * t643;
t610 = t539 * t643;
t609 = t537 * t641;
t608 = t540 * t641;
t607 = t538 * t639;
t606 = t541 * t639;
t605 = t473 * t635;
t604 = t474 * t632;
t603 = t475 * t629;
t503 = t506 ^ 2;
t518 = pkin(1) * t555 + pkin(2) * t533;
t602 = t503 * t518 / 0.4e1;
t601 = t503 * t658;
t504 = t507 ^ 2;
t519 = pkin(1) * t557 + pkin(2) * t534;
t600 = t504 * t519 / 0.4e1;
t599 = t504 * t658;
t505 = t508 ^ 2;
t520 = pkin(1) * t559 + pkin(2) * t535;
t598 = t505 * t520 / 0.4e1;
t597 = t505 * t658;
t593 = t647 * t671;
t592 = t646 * t670;
t591 = t645 * t669;
t590 = t515 * t617;
t589 = t516 * t616;
t588 = t517 * t615;
t584 = -t515 * t571 * t593 + (t601 * t665 + t473) * t509;
t583 = -t516 * t571 * t592 + (t599 * t664 + t474) * t510;
t582 = -t517 * t571 * t591 + (t597 * t663 + t475) * t511;
t581 = -t515 * t491 * t614 + (t543 * t601 - t605) * t509;
t580 = -t516 * t492 * t613 + (t546 * t599 - t604) * t510;
t579 = -t517 * t493 * t612 + (t549 * t597 - t603) * t511;
t1 = [0, (t473 * t611 + t474 * t609 + t475 * t607) * t571, 0, 0, (t467 * t611 + t468 * t609 + t469 * t607 + (-t536 * t620 - t537 * t619 - t538 * t618) * t570) * t571, t536 * t590 + t537 * t589 + t538 * t588 + (t581 * t536 + t580 * t537 + t579 * t538) * t570, -t536 * t644 - t537 * t642 - t538 * t640 + (t584 * t536 + t583 * t537 + t582 * t538) * t570, 0; 0, (t473 * t610 + t474 * t608 + t475 * t606) * t571, 0, 0, (t467 * t610 + t468 * t608 + t469 * t606 + (-t539 * t620 - t540 * t619 - t541 * t618) * t570) * t571, t539 * t590 + t540 * t589 + t541 * t588 + (t581 * t539 + t580 * t540 + t579 * t541) * t570, -t539 * t644 - t540 * t642 - t541 * t640 + (t584 * t539 + t583 * t540 + t582 * t541) * t570, 0; 0, (-t473 * t638 - t474 * t637 - t475 * t636) * t571, 0, 0, (-t467 * t638 - t468 * t637 - t469 * t636 + (t518 * t653 + t519 * t652 + t520 * t651) * t570) * t571, -t533 * t617 - t534 * t616 - t535 * t615 + (t518 * t605 + t519 * t604 + t520 * t603 + (t543 * t602 + t546 * t600 + t549 * t598) * t572 + (t636 * t645 + t637 * t646 + t638 * t647) * t571) * t570, t533 * t470 + t534 * t471 + t535 * t472 + (-t518 * t473 - t519 * t474 - t520 * t475 + (t598 * t663 + t600 * t664 + t602 * t665) * t572 + (t533 * t593 + t534 * t592 + t535 * t591) * t571) * t570, 0;];
tau_reg  = t1;
