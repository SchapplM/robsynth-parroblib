% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRR2G3P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR2G3P1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR2G3P1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3P1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G3P1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR2G3P1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3P1A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G3P1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3P1A0_invdyn_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3P1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3P1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR2G3P1A0_invdyn_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:07
% EndTime: 2020-03-09 21:20:08
% DurationCPUTime: 1.06s
% Computational Cost: add. (5045->172), mult. (3930->310), div. (1710->9), fcn. (2844->30), ass. (0->171)
t723 = -2 * pkin(1);
t722 = 2 * pkin(1);
t721 = 2 * pkin(2);
t649 = cos(qJ(3,3));
t720 = t649 * pkin(1);
t651 = cos(qJ(3,2));
t719 = t651 * pkin(1);
t653 = cos(qJ(3,1));
t718 = t653 * pkin(1);
t638 = legFrame(3,2);
t622 = -t638 + qJ(2,3);
t616 = qJ(3,3) + t622;
t604 = sin(t616);
t607 = cos(t616);
t655 = xDP(2);
t656 = xDP(1);
t663 = t604 * t656 - t607 * t655;
t643 = sin(qJ(3,3));
t632 = 0.1e1 / t643;
t659 = 1 / pkin(1);
t687 = t632 * t659;
t571 = t663 * t687;
t610 = sin(t622);
t613 = cos(t622);
t562 = t663 * pkin(2) + (t610 * t656 - t613 * t655) * pkin(1);
t658 = 1 / pkin(2);
t679 = t658 * t659;
t672 = t562 * t679;
t666 = t632 * t672;
t559 = -t571 + t666;
t657 = pkin(2) ^ 2;
t708 = t663 * t632;
t556 = -t571 + t666 / 0.2e1;
t714 = t556 * t649;
t717 = (t559 * t657 + (t714 * t721 - t708) * pkin(1)) * t663;
t639 = legFrame(2,2);
t623 = -t639 + qJ(2,2);
t617 = qJ(3,2) + t623;
t605 = sin(t617);
t608 = cos(t617);
t662 = t605 * t656 - t608 * t655;
t645 = sin(qJ(3,2));
t634 = 0.1e1 / t645;
t685 = t634 * t659;
t572 = t662 * t685;
t611 = sin(t623);
t614 = cos(t623);
t563 = t662 * pkin(2) + (t611 * t656 - t614 * t655) * pkin(1);
t671 = t563 * t679;
t665 = t634 * t671;
t560 = -t572 + t665;
t707 = t662 * t634;
t557 = -t572 + t665 / 0.2e1;
t713 = t557 * t651;
t716 = (t560 * t657 + (t713 * t721 - t707) * pkin(1)) * t662;
t640 = legFrame(1,2);
t624 = -t640 + qJ(2,1);
t618 = qJ(3,1) + t624;
t606 = sin(t618);
t609 = cos(t618);
t661 = t606 * t656 - t609 * t655;
t647 = sin(qJ(3,1));
t636 = 0.1e1 / t647;
t683 = t636 * t659;
t573 = t661 * t683;
t612 = sin(t624);
t615 = cos(t624);
t564 = t661 * pkin(2) + (t612 * t656 - t615 * t655) * pkin(1);
t670 = t564 * t679;
t664 = t636 * t670;
t561 = -t573 + t664;
t706 = t661 * t636;
t558 = -t573 + t664 / 0.2e1;
t712 = t558 * t653;
t715 = (t561 * t657 + (t712 * t721 - t706) * pkin(1)) * t661;
t711 = t559 * t562;
t710 = t560 * t563;
t709 = t561 * t564;
t580 = pkin(1) * t610 + pkin(2) * t604;
t705 = t580 * t632;
t642 = xDDP(1);
t704 = t580 * t642;
t581 = pkin(1) * t611 + pkin(2) * t605;
t703 = t581 * t634;
t702 = t581 * t642;
t582 = pkin(1) * t612 + pkin(2) * t606;
t701 = t582 * t636;
t700 = t582 * t642;
t583 = -pkin(1) * t613 - pkin(2) * t607;
t699 = t583 * t632;
t641 = xDDP(2);
t698 = t583 * t641;
t584 = -pkin(1) * t614 - pkin(2) * t608;
t697 = t584 * t634;
t696 = t584 * t641;
t585 = -pkin(1) * t615 - pkin(2) * t609;
t695 = t585 * t636;
t694 = t585 * t641;
t693 = t604 * t632;
t692 = t605 * t634;
t691 = t606 * t636;
t690 = t607 * t632;
t689 = t608 * t634;
t688 = t609 * t636;
t633 = 0.1e1 / t643 ^ 2;
t660 = 1 / pkin(1) ^ 2;
t686 = t633 * t660;
t635 = 0.1e1 / t645 ^ 2;
t684 = t635 * t660;
t637 = 0.1e1 / t647 ^ 2;
t682 = t637 * t660;
t681 = t641 * t659;
t680 = t642 * t659;
t678 = g(1) * t607 + g(2) * t604;
t677 = g(1) * t608 + g(2) * t605;
t676 = g(1) * t609 + g(2) * t606;
t675 = (pkin(2) + t720) * t711;
t674 = (pkin(2) + t719) * t710;
t673 = (pkin(2) + t718) * t709;
t669 = -g(1) * t604 + g(2) * t607;
t668 = -g(1) * t605 + g(2) * t608;
t667 = -g(1) * t606 + g(2) * t609;
t544 = -t680 * t693 + t681 * t690 + ((-pkin(2) * t559 + t649 * t708) * t663 + t711) * t686;
t545 = -t680 * t692 + t681 * t689 + ((-pkin(2) * t560 + t651 * t707) * t662 + t710) * t684;
t546 = -t680 * t691 + t681 * t688 + ((-t561 * pkin(2) + t653 * t706) * t661 + t709) * t682;
t654 = cos(qJ(2,1));
t652 = cos(qJ(2,2));
t650 = cos(qJ(2,3));
t648 = sin(qJ(2,1));
t646 = sin(qJ(2,2));
t644 = sin(qJ(2,3));
t631 = cos(t640);
t630 = cos(t639);
t629 = cos(t638);
t628 = sin(t640);
t627 = sin(t639);
t626 = sin(t638);
t594 = t631 * g(1) - t628 * g(2);
t593 = t630 * g(1) - t627 * g(2);
t592 = t629 * g(1) - t626 * g(2);
t591 = t628 * g(1) + t631 * g(2);
t590 = t627 * g(1) + t630 * g(2);
t589 = t626 * g(1) + t629 * g(2);
t576 = t661 ^ 2;
t575 = t662 ^ 2;
t574 = t663 ^ 2;
t570 = t591 * t654 - t594 * t648;
t569 = t591 * t648 + t594 * t654;
t568 = t590 * t652 - t593 * t646;
t567 = t590 * t646 + t593 * t652;
t566 = t589 * t650 - t592 * t644;
t565 = t589 * t644 + t592 * t650;
t543 = t546 * t718 + t576 * t683 + t676;
t542 = t544 * t720 + t574 * t687 + t678;
t541 = t545 * t719 + t575 * t685 + t677;
t540 = t653 * t576 * t637 * t659 - t647 * t546 * pkin(1) + t667;
t539 = t651 * t575 * t635 * t659 - t645 * t545 * pkin(1) + t668;
t538 = t649 * t574 * t659 * t633 - t643 * t544 * pkin(1) + t669;
t537 = ((-t673 + t715) * t682 + (t694 + t700) * t683) * t658 + t546;
t536 = ((-t674 + t716) * t684 + (t696 + t702) * t685) * t658 + t545;
t535 = ((-t675 + t717) * t686 + (t698 + t704) * t687) * t658 + t544;
t534 = ((t715 / 0.2e1 - t673 / 0.2e1) * t682 + (t700 / 0.2e1 + t694 / 0.2e1) * t683) * t658 + t546;
t533 = ((t716 / 0.2e1 - t674 / 0.2e1) * t684 + (t702 / 0.2e1 + t696 / 0.2e1) * t685) * t658 + t545;
t532 = ((t717 / 0.2e1 - t675 / 0.2e1) * t686 + (t704 / 0.2e1 + t698 / 0.2e1) * t687) * t658 + t544;
t531 = (t534 * t653 - t558 * t670) * t722 + t676;
t530 = (t533 * t651 - t557 * t671) * t722 + t677;
t529 = (t532 * t649 - t556 * t672) * t722 + t678;
t528 = (t647 * t534 + t664 * t712) * t723 + t667;
t527 = (t645 * t533 + t665 * t713) * t723 + t668;
t526 = (t643 * t532 + t666 * t714) * t723 + t669;
t1 = [(t642 - g(1)) * MDP(8) + ((-t544 * t693 - t545 * t692 - t546 * t691) * MDP(2) + (-t565 * t693 - t567 * t692 - t569 * t691) * MDP(3) + (-t566 * t693 - t568 * t692 - t570 * t691) * MDP(4) + (-t535 * t693 - t536 * t692 - t537 * t691) * MDP(5) + (-t529 * t693 - t530 * t692 - t531 * t691) * MDP(6) + (-t526 * t693 - t527 * t692 - t528 * t691) * MDP(7) + ((t535 * t705 + t536 * t703 + t537 * t701) * MDP(5) + (t541 * t703 + t542 * t705 + t543 * t701) * MDP(6) + (t538 * t705 + t539 * t703 + t540 * t701) * MDP(7)) * t658) * t659; (t641 - g(2)) * MDP(8) + ((t544 * t690 + t545 * t689 + t546 * t688) * MDP(2) + (t565 * t690 + t567 * t689 + t569 * t688) * MDP(3) + (t566 * t690 + t568 * t689 + t570 * t688) * MDP(4) + (t535 * t690 + t536 * t689 + t537 * t688) * MDP(5) + (t529 * t690 + t530 * t689 + t531 * t688) * MDP(6) + (t526 * t690 + t527 * t689 + t528 * t688) * MDP(7) + ((t535 * t699 + t536 * t697 + t537 * t695) * MDP(5) + (t541 * t697 + t542 * t699 + t543 * t695) * MDP(6) + (t538 * t699 + t539 * t697 + t540 * t695) * MDP(7)) * t658) * t659; ((3 * MDP(1)) + MDP(8)) * (xDDP(3) - g(3));];
tauX  = t1;
