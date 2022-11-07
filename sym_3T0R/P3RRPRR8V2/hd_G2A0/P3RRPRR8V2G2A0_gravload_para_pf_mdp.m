% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V2G2A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V2G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V2G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:08:20
% EndTime: 2022-11-07 13:08:22
% DurationCPUTime: 2.33s
% Computational Cost: add. (1227->227), mult. (2067->376), div. (165->9), fcn. (1815->26), ass. (0->183)
t4694 = qJ(3,1) + pkin(5);
t4682 = pkin(6) + t4694;
t4666 = 0.1e1 / t4682;
t4709 = cos(qJ(1,1));
t4780 = t4666 * t4709;
t4703 = sin(qJ(1,1));
t4697 = legFrame(1,2);
t4669 = sin(t4697);
t4672 = cos(t4697);
t4721 = t4672 * g(1) - t4669 * g(2);
t4810 = g(3) * t4703 - t4721 * t4709;
t4819 = t4810 * t4780;
t4691 = cos(pkin(7));
t4795 = t4691 * pkin(3);
t4653 = pkin(2) + t4795;
t4708 = cos(qJ(2,1));
t4690 = sin(pkin(7));
t4702 = sin(qJ(2,1));
t4768 = t4690 * t4702;
t4757 = pkin(3) * t4768;
t4792 = 0.1e1 / (t4653 * t4708 - t4757) * t4666;
t4818 = t4810 * t4792;
t4693 = qJ(3,2) + pkin(5);
t4681 = pkin(6) + t4693;
t4665 = 0.1e1 / t4681;
t4707 = cos(qJ(1,2));
t4781 = t4665 * t4707;
t4701 = sin(qJ(1,2));
t4696 = legFrame(2,2);
t4668 = sin(t4696);
t4671 = cos(t4696);
t4722 = t4671 * g(1) - t4668 * g(2);
t4811 = g(3) * t4701 - t4722 * t4707;
t4817 = t4811 * t4781;
t4706 = cos(qJ(2,2));
t4700 = sin(qJ(2,2));
t4769 = t4690 * t4700;
t4758 = pkin(3) * t4769;
t4793 = 0.1e1 / (t4653 * t4706 - t4758) * t4665;
t4816 = t4811 * t4793;
t4692 = qJ(3,3) + pkin(5);
t4680 = pkin(6) + t4692;
t4664 = 0.1e1 / t4680;
t4705 = cos(qJ(1,3));
t4782 = t4664 * t4705;
t4699 = sin(qJ(1,3));
t4695 = legFrame(3,2);
t4667 = sin(t4695);
t4670 = cos(t4695);
t4723 = t4670 * g(1) - t4667 * g(2);
t4812 = g(3) * t4699 - t4723 * t4705;
t4815 = t4812 * t4782;
t4704 = cos(qJ(2,3));
t4698 = sin(qJ(2,3));
t4770 = t4690 * t4698;
t4759 = pkin(3) * t4770;
t4794 = 0.1e1 / (t4653 * t4704 - t4759) * t4664;
t4814 = t4812 * t4794;
t4813 = MDP(3) - MDP(13);
t4809 = 0.2e1 * t4704 ^ 2;
t4808 = 0.2e1 * t4706 ^ 2;
t4807 = 0.2e1 * t4708 ^ 2;
t4806 = MDP(14) * pkin(2);
t4684 = qJ(2,3) + pkin(7);
t4661 = cos(t4684);
t4805 = pkin(3) * t4661;
t4685 = qJ(2,2) + pkin(7);
t4662 = cos(t4685);
t4804 = pkin(3) * t4662;
t4686 = qJ(2,1) + pkin(7);
t4663 = cos(t4686);
t4803 = pkin(3) * t4663;
t4738 = pkin(1) * t4699 - t4705 * t4680;
t4683 = t4691 ^ 2;
t4766 = pkin(3) * (t4683 - 0.1e1);
t4799 = (t4699 * t4766 + t4738 * t4770) * pkin(3);
t4737 = pkin(1) * t4701 - t4707 * t4681;
t4798 = (t4701 * t4766 + t4737 * t4769) * pkin(3);
t4736 = pkin(1) * t4703 - t4709 * t4682;
t4797 = (t4703 * t4766 + t4736 * t4768) * pkin(3);
t4796 = t4690 * pkin(3);
t4673 = t4704 * pkin(2);
t4644 = 0.1e1 / (t4673 + t4805);
t4791 = t4644 * t4667;
t4790 = t4644 * t4670;
t4674 = t4706 * pkin(2);
t4645 = 0.1e1 / (t4674 + t4804);
t4789 = t4645 * t4668;
t4788 = t4645 * t4671;
t4675 = t4708 * pkin(2);
t4646 = 0.1e1 / (t4675 + t4803);
t4787 = t4646 * t4669;
t4786 = t4646 * t4672;
t4785 = t4653 * t4670;
t4784 = t4653 * t4671;
t4783 = t4653 * t4672;
t4779 = t4667 * t4653;
t4778 = t4667 * t4699;
t4777 = t4668 * t4653;
t4776 = t4668 * t4701;
t4775 = t4669 * t4653;
t4774 = t4669 * t4703;
t4773 = t4670 * t4699;
t4772 = t4671 * t4701;
t4771 = t4672 * t4703;
t4767 = pkin(2) * t4795;
t4765 = t4670 * t4796;
t4764 = t4671 * t4796;
t4763 = t4672 * t4796;
t4762 = t4667 * t4796;
t4761 = t4668 * t4796;
t4760 = t4669 * t4796;
t4719 = g(3) * t4705 + t4723 * t4699;
t4747 = t4719 * t4794;
t4718 = g(3) * t4707 + t4722 * t4701;
t4746 = t4718 * t4793;
t4717 = g(3) * t4709 + t4721 * t4703;
t4745 = t4717 * t4792;
t4654 = t4673 + pkin(1);
t4651 = t4654 * t4705;
t4598 = (t4699 * t4654 - t4692 * t4705) * g(3) - t4723 * (t4699 * t4692 + t4651);
t4744 = t4598 * t4794;
t4655 = t4674 + pkin(1);
t4652 = t4655 * t4707;
t4599 = (t4701 * t4655 - t4693 * t4707) * g(3) - (t4701 * t4693 + t4652) * t4722;
t4743 = t4599 * t4793;
t4656 = t4675 + pkin(1);
t4650 = t4709 * t4656;
t4600 = (t4703 * t4656 - t4694 * t4709) * g(3) - (t4703 * t4694 + t4650) * t4721;
t4742 = t4600 * t4792;
t4735 = t4661 * t4814;
t4734 = t4704 * t4814;
t4733 = t4662 * t4816;
t4732 = t4706 * t4816;
t4731 = t4663 * t4818;
t4730 = t4708 * t4818;
t4658 = sin(t4684);
t4729 = t4658 * t4814;
t4728 = t4698 * t4814;
t4659 = sin(t4685);
t4727 = t4659 * t4816;
t4726 = t4700 * t4816;
t4660 = sin(t4686);
t4725 = t4660 * t4818;
t4724 = t4702 * t4818;
t4710 = pkin(3) ^ 2;
t4711 = pkin(2) ^ 2;
t4720 = 0.2e1 * t4683 * t4710 - t4710 + t4711 + 0.2e1 * t4767;
t4641 = t4667 * g(1) + t4670 * g(2);
t4607 = -t4641 * t4704 + t4698 * t4719;
t4642 = t4668 * g(1) + t4671 * g(2);
t4609 = -t4642 * t4706 + t4700 * t4718;
t4643 = t4669 * g(1) + t4672 * g(2);
t4611 = -t4643 * t4708 + t4702 * t4717;
t4716 = t4607 * t4791 + t4609 * t4789 + t4611 * t4787;
t4715 = t4607 * t4790 + t4609 * t4788 + t4611 * t4786;
t4657 = pkin(1) * t4796;
t4649 = pkin(1) * t4702 - t4796;
t4648 = pkin(1) * t4700 - t4796;
t4647 = pkin(1) * t4698 - t4796;
t4640 = t4767 + t4711 / 0.2e1 + (t4683 - 0.1e1 / 0.2e1) * t4710;
t4621 = -0.2e1 * t4703 * t4757 + t4736;
t4620 = -0.2e1 * t4701 * t4758 + t4737;
t4619 = -0.2e1 * t4699 * t4759 + t4738;
t4618 = t4720 * t4702 + t4657;
t4617 = t4720 * t4700 + t4657;
t4616 = t4720 * t4698 + t4657;
t4612 = t4643 * t4702 + t4708 * t4717;
t4610 = t4642 * t4700 + t4706 * t4718;
t4608 = t4641 * t4698 + t4704 * t4719;
t4606 = t4643 * t4660 + t4663 * t4717;
t4605 = -t4643 * t4663 + t4660 * t4717;
t4604 = t4642 * t4659 + t4662 * t4718;
t4603 = -t4642 * t4662 + t4659 * t4718;
t4602 = t4641 * t4658 + t4661 * t4719;
t4601 = -t4641 * t4661 + t4658 * t4719;
t4597 = (t4653 * t4771 + t4760) * t4708 + t4702 * (-t4703 * t4763 + t4775);
t4596 = (t4653 * t4772 + t4761) * t4706 + t4700 * (-t4701 * t4764 + t4777);
t4595 = (t4653 * t4773 + t4762) * t4704 + t4698 * (-t4699 * t4765 + t4779);
t4594 = (-t4653 * t4774 + t4763) * t4708 + (t4703 * t4760 + t4783) * t4702;
t4593 = (-t4653 * t4776 + t4764) * t4706 + (t4701 * t4761 + t4784) * t4700;
t4592 = (-t4653 * t4778 + t4765) * t4704 + (t4699 * t4762 + t4785) * t4698;
t1 = [(t4595 * t4814 + t4596 * t4816 + t4597 * t4818) * MDP(2) + (t4595 * t4734 + t4596 * t4732 + t4597 * t4730 + t4716) * MDP(9) + (-t4595 * t4728 - t4596 * t4726 - t4597 * t4724 + t4608 * t4791 + t4610 * t4789 + t4612 * t4787) * MDP(10) + (t4595 * t4735 + t4596 * t4733 + t4597 * t4731 + t4601 * t4791 + t4603 * t4789 + t4605 * t4787) * MDP(11) + (-t4595 * t4729 - t4596 * t4727 - t4597 * t4725 + t4602 * t4791 + t4604 * t4789 + t4606 * t4787) * MDP(12) + (t4597 * t4742 - ((t4640 * t4771 + t4653 * t4760) * t4807 + (t4618 * t4669 + t4621 * t4783) * t4708 - t4672 * t4797 + t4649 * t4775) * t4818 + t4596 * t4743 - ((t4640 * t4772 + t4653 * t4761) * t4808 + (t4617 * t4668 + t4620 * t4784) * t4706 - t4671 * t4798 + t4648 * t4777) * t4816 + t4595 * t4744 - ((t4640 * t4773 + t4653 * t4762) * t4809 + (t4616 * t4667 + t4619 * t4785) * t4704 - t4670 * t4799 + t4647 * t4779) * t4814) * MDP(14) - g(1) * MDP(15) + t4716 * t4806 + t4813 * (t4595 * t4747 + t4596 * t4746 + t4597 * t4745); (t4592 * t4814 + t4593 * t4816 + t4594 * t4818) * MDP(2) + (t4592 * t4734 + t4593 * t4732 + t4594 * t4730 + t4715) * MDP(9) + (-t4592 * t4728 - t4593 * t4726 - t4594 * t4724 + t4608 * t4790 + t4610 * t4788 + t4612 * t4786) * MDP(10) + (t4592 * t4735 + t4593 * t4733 + t4594 * t4731 + t4601 * t4790 + t4603 * t4788 + t4605 * t4786) * MDP(11) + (-t4592 * t4729 - t4593 * t4727 - t4594 * t4725 + t4602 * t4790 + t4604 * t4788 + t4606 * t4786) * MDP(12) + (t4594 * t4742 - ((-t4640 * t4774 + t4653 * t4763) * t4807 + (t4672 * t4618 - t4621 * t4775) * t4708 + t4669 * t4797 + t4649 * t4783) * t4818 + t4593 * t4743 - ((-t4640 * t4776 + t4653 * t4764) * t4808 + (t4671 * t4617 - t4620 * t4777) * t4706 + t4668 * t4798 + t4648 * t4784) * t4816 + t4592 * t4744 - ((-t4640 * t4778 + t4653 * t4765) * t4809 + (t4670 * t4616 - t4619 * t4779) * t4704 + t4667 * t4799 + t4647 * t4785) * t4814) * MDP(14) - g(2) * MDP(15) + t4715 * t4806 + t4813 * (t4592 * t4747 + t4593 * t4746 + t4594 * t4745); (t4819 + t4817 + t4815) * MDP(2) + (t4704 * t4815 + t4706 * t4817 + t4708 * t4819) * MDP(9) + (-t4698 * t4815 - t4700 * t4817 - t4702 * t4819) * MDP(10) + (t4661 * t4815 + t4662 * t4817 + t4663 * t4819) * MDP(11) + (-t4658 * t4815 - t4659 * t4817 - t4660 * t4819) * MDP(12) + ((t4709 * t4600 - (t4703 * t4682 + t4709 * t4803 + t4650) * t4810) * t4666 + (t4707 * t4599 - (t4701 * t4681 + t4707 * t4804 + t4652) * t4811) * t4665 + (t4705 * t4598 - (t4699 * t4680 + t4705 * t4805 + t4651) * t4812) * t4664) * MDP(14) - g(3) * MDP(15) + t4813 * (t4717 * t4780 + t4718 * t4781 + t4719 * t4782);];
taugX  = t1;
