% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRR1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR1G2A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G2A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2A0_coriolisvec_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:28
% EndTime: 2020-03-09 21:18:30
% DurationCPUTime: 1.95s
% Computational Cost: add. (19368->159), mult. (11652->358), div. (2469->9), fcn. (11733->24), ass. (0->177)
t610 = legFrame(3,2);
t601 = sin(t610);
t604 = cos(t610);
t620 = xDP(2);
t621 = xDP(1);
t574 = -t601 * t620 + t604 * t621;
t607 = pkin(7) + qJ(2,3);
t598 = qJ(3,3) + t607;
t586 = sin(t598);
t589 = cos(t598);
t619 = xDP(3);
t637 = -t574 * t586 - t589 * t619;
t624 = 0.1e1 / pkin(2);
t592 = sin(t607);
t595 = cos(t607);
t568 = -t595 * t586 + t589 * t592;
t722 = 0.1e1 / t568;
t687 = t722 * t624;
t544 = -pkin(2) * (t574 * t592 + t595 * t619) + t637 * pkin(3);
t623 = 0.1e1 / pkin(3);
t704 = t544 * t623;
t539 = (-t637 + t704) * t687;
t634 = t586 * t592 + t589 * t595;
t698 = t637 * t722;
t532 = -t539 * pkin(3) + t634 * t698;
t625 = 0.1e1 / pkin(2) ^ 2;
t560 = 0.1e1 / t568 ^ 2;
t644 = t539 * (t634 * pkin(2) + pkin(3)) * t560 * t704;
t695 = t722 ^ 2;
t709 = t539 * t544;
t622 = pkin(3) ^ 2;
t662 = 0.2e1 * pkin(3) * t624;
t713 = (t539 * t622 + (-t698 + t634 * (-t637 + t704 / 0.2e1) * t722 * t662) * pkin(2)) * t623;
t520 = (t644 - (t709 - (-t532 - t713) * t637) * t695) * t625;
t577 = pkin(2) * t592 + pkin(3) * t586;
t694 = t722 * t577;
t728 = t520 * t694;
t611 = legFrame(2,2);
t602 = sin(t611);
t605 = cos(t611);
t575 = -t602 * t620 + t605 * t621;
t608 = pkin(7) + qJ(2,2);
t599 = qJ(3,2) + t608;
t587 = sin(t599);
t590 = cos(t599);
t636 = -t575 * t587 - t590 * t619;
t593 = sin(t608);
t596 = cos(t608);
t569 = -t596 * t587 + t590 * t593;
t721 = 0.1e1 / t569;
t684 = t721 * t624;
t545 = -pkin(2) * (t575 * t593 + t596 * t619) + t636 * pkin(3);
t703 = t545 * t623;
t541 = (-t636 + t703) * t684;
t633 = t587 * t593 + t590 * t596;
t697 = t636 * t721;
t533 = -t541 * pkin(3) + t633 * t697;
t563 = 0.1e1 / t569 ^ 2;
t641 = t541 * (t633 * pkin(2) + pkin(3)) * t563 * t703;
t693 = t721 ^ 2;
t707 = t541 * t545;
t712 = (t541 * t622 + (-t697 + t633 * (-t636 + t703 / 0.2e1) * t721 * t662) * pkin(2)) * t623;
t522 = (t641 - (t707 - (-t533 - t712) * t636) * t693) * t625;
t578 = pkin(2) * t593 + pkin(3) * t587;
t692 = t721 * t578;
t727 = t522 * t692;
t612 = legFrame(1,2);
t603 = sin(t612);
t606 = cos(t612);
t576 = -t603 * t620 + t606 * t621;
t609 = pkin(7) + qJ(2,1);
t600 = qJ(3,1) + t609;
t588 = sin(t600);
t591 = cos(t600);
t635 = -t576 * t588 - t591 * t619;
t594 = sin(t609);
t597 = cos(t609);
t570 = -t597 * t588 + t591 * t594;
t720 = 0.1e1 / t570;
t681 = t720 * t624;
t546 = -pkin(2) * (t576 * t594 + t597 * t619) + t635 * pkin(3);
t702 = t546 * t623;
t543 = (-t635 + t702) * t681;
t632 = t588 * t594 + t591 * t597;
t696 = t635 * t720;
t534 = -t543 * pkin(3) + t632 * t696;
t566 = 0.1e1 / t570 ^ 2;
t638 = t543 * (t632 * pkin(2) + pkin(3)) * t566 * t702;
t691 = t720 ^ 2;
t705 = t543 * t546;
t711 = (t543 * t622 + (-t696 + t632 * (-t635 + t702 / 0.2e1) * t720 * t662) * pkin(2)) * t623;
t524 = (t638 - (t705 - (-t534 - t711) * t635) * t691) * t625;
t579 = pkin(2) * t594 + pkin(3) * t588;
t690 = t720 * t579;
t726 = t524 * t690;
t725 = t586 * t722;
t724 = t587 * t721;
t723 = t588 * t720;
t526 = (t532 * t637 + t709) * t625 * t695;
t613 = sin(qJ(3,3));
t719 = t526 * t613;
t616 = cos(qJ(3,3));
t718 = t526 * t616;
t527 = (t533 * t636 + t707) * t625 * t693;
t614 = sin(qJ(3,2));
t717 = t527 * t614;
t617 = cos(qJ(3,2));
t716 = t527 * t617;
t528 = (t534 * t635 + t705) * t625 * t691;
t615 = sin(qJ(3,1));
t715 = t528 * t615;
t618 = cos(qJ(3,1));
t714 = t528 * t618;
t710 = (-0.2e1 * t637 + t704) * t687 * t544;
t708 = (-0.2e1 * t636 + t703) * t684 * t545;
t706 = (-0.2e1 * t635 + t702) * t681 * t546;
t547 = t637 ^ 2;
t701 = t547 * t625;
t548 = t636 ^ 2;
t700 = t548 * t625;
t549 = t635 ^ 2;
t699 = t549 * t625;
t580 = -pkin(2) * t595 - pkin(3) * t589;
t689 = t722 * t580;
t688 = t722 * t589;
t581 = -pkin(2) * t596 - pkin(3) * t590;
t686 = t721 * t581;
t685 = t721 * t590;
t582 = -pkin(2) * t597 - pkin(3) * t591;
t683 = t720 * t582;
t682 = t720 * t591;
t680 = t526 * t689;
t679 = t527 * t686;
t678 = t528 * t683;
t677 = t560 * t701;
t561 = t722 * t560;
t676 = t547 * t561 * t580;
t675 = t563 * t700;
t564 = t721 * t563;
t674 = t548 * t564 * t581;
t673 = t566 * t699;
t567 = t720 * t566;
t672 = t549 * t567 * t582;
t671 = t601 * t725;
t670 = t602 * t724;
t669 = t603 * t723;
t668 = t604 * t725;
t667 = t605 * t724;
t666 = t606 * t723;
t521 = (t644 - (0.2e1 * t709 - (-0.2e1 * t532 - t713) * t637) * t695) * t625;
t665 = t521 * t688;
t523 = (t641 - (0.2e1 * t707 - (-0.2e1 * t533 - t712) * t636) * t693) * t625;
t664 = t523 * t685;
t525 = (t638 - (0.2e1 * t705 - (-0.2e1 * t534 - t711) * t635) * t691) * t625;
t663 = t525 * t682;
t661 = t560 * t589 * t710;
t660 = t586 * t624 * t710;
t659 = t563 * t590 * t708;
t658 = t587 * t624 * t708;
t657 = t566 * t591 * t706;
t656 = t588 * t624 * t706;
t655 = t561 * t577 * t701;
t654 = t564 * t578 * t700;
t653 = t567 * t579 * t699;
t652 = t521 * t671;
t651 = t523 * t670;
t650 = t525 * t669;
t649 = t521 * t668;
t648 = t523 * t667;
t647 = t525 * t666;
t646 = t613 * t660;
t645 = t616 * t660;
t643 = t614 * t658;
t642 = t617 * t658;
t640 = t615 * t656;
t639 = t618 * t656;
t1 = [0, (t526 * t668 + t527 * t667 + t528 * t666) * t624, 0, 0, (-t520 * t668 - t522 * t667 - t524 * t666 + (t604 * t728 + t605 * t727 + t606 * t726) * t623) * t624, -t616 * t649 - t617 * t648 - t618 * t647 + ((t566 * t640 + (t615 * t673 - t714) * t690) * t606 + (t563 * t643 + (t614 * t675 - t716) * t692) * t605 + (t560 * t646 + (t613 * t677 - t718) * t694) * t604) * t623, t613 * t649 + t614 * t648 + t615 * t647 + ((t566 * t639 + (t618 * t673 + t715) * t690) * t606 + (t563 * t642 + (t617 * t675 + t717) * t692) * t605 + (t560 * t645 + (t616 * t677 + t719) * t694) * t604) * t623, 0; 0, (-t526 * t671 - t527 * t670 - t528 * t669) * t624, 0, 0, (t520 * t671 + t522 * t670 + t524 * t669 + (-t601 * t728 - t602 * t727 - t603 * t726) * t623) * t624, t616 * t652 + t617 * t651 + t618 * t650 + ((-t615 * t653 - (-t579 * t714 + t640 * t720) * t720) * t603 + (-t614 * t654 - (-t578 * t716 + t643 * t721) * t721) * t602 + (-t613 * t655 - (-t577 * t718 + t646 * t722) * t722) * t601) * t623, -t613 * t652 - t614 * t651 - t615 * t650 + ((-t618 * t653 - (t579 * t715 + t639 * t720) * t720) * t603 + (-t617 * t654 - (t578 * t717 + t642 * t721) * t721) * t602 + (-t616 * t655 - (t577 * t719 + t645 * t722) * t722) * t601) * t623, 0; 0, (t526 * t688 + t527 * t685 + t528 * t682) * t624, 0, 0, (-t520 * t688 - t522 * t685 - t524 * t682 + (-t520 * t689 - t522 * t686 - t524 * t683) * t623) * t624, -t616 * t665 - t617 * t664 - t618 * t663 + (t616 * t680 + t617 * t679 + t618 * t678 + (-t613 * t676 - t614 * t674 - t615 * t672) * t625 + (t613 * t661 + t614 * t659 + t615 * t657) * t624) * t623, t613 * t665 + t614 * t664 + t615 * t663 + (-t613 * t680 - t614 * t679 - t615 * t678 + (-t616 * t676 - t617 * t674 - t618 * t672) * t625 + (t616 * t661 + t617 * t659 + t618 * t657) * t624) * t623, 0;];
tau_reg  = t1;
