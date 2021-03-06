% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRRRR7V1G3A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [18x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRRRR7V1G3A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 09:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RRRRR7V1G3A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(18,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_mdp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 09:00:29
% EndTime: 2020-08-07 09:00:50
% DurationCPUTime: 22.04s
% Computational Cost: add. (18558->799), mult. (32379->1360), div. (3480->27), fcn. (26529->60), ass. (0->581)
t4576 = sin(qJ(1,3));
t4592 = -pkin(5) - pkin(4);
t4532 = t4576 * t4592;
t4585 = cos(qJ(1,3));
t4583 = cos(qJ(3,3));
t4940 = pkin(2) * t4583;
t4527 = pkin(1) + t4940;
t4584 = cos(qJ(2,3));
t4574 = sin(qJ(3,3));
t4575 = sin(qJ(2,3));
t4870 = t4574 * t4575;
t4655 = pkin(2) * t4870 - t4527 * t4584;
t4953 = t4655 * t4585 + t4532;
t4579 = sin(qJ(1,2));
t4533 = t4579 * t4592;
t4588 = cos(qJ(1,2));
t4586 = cos(qJ(3,2));
t4939 = pkin(2) * t4586;
t4529 = pkin(1) + t4939;
t4587 = cos(qJ(2,2));
t4577 = sin(qJ(3,2));
t4578 = sin(qJ(2,2));
t4866 = t4577 * t4578;
t4654 = pkin(2) * t4866 - t4529 * t4587;
t4952 = t4654 * t4588 + t4533;
t4582 = sin(qJ(1,1));
t4534 = t4582 * t4592;
t4591 = cos(qJ(1,1));
t4589 = cos(qJ(3,1));
t4938 = pkin(2) * t4589;
t4531 = pkin(1) + t4938;
t4590 = cos(qJ(2,1));
t4580 = sin(qJ(3,1));
t4581 = sin(qJ(2,1));
t4862 = t4580 * t4581;
t4653 = pkin(2) * t4862 - t4531 * t4590;
t4951 = t4653 * t4591 + t4534;
t4571 = legFrame(3,2);
t4541 = sin(t4571);
t4544 = cos(t4571);
t4950 = t4541 * t4544;
t4572 = legFrame(2,2);
t4542 = sin(t4572);
t4545 = cos(t4572);
t4949 = t4542 * t4545;
t4573 = legFrame(1,2);
t4543 = sin(t4573);
t4546 = cos(t4573);
t4948 = t4543 * t4546;
t4947 = -0.2e1 * pkin(1);
t4946 = 0.2e1 * pkin(1);
t4945 = 0.2e1 * pkin(2);
t4944 = 2 * MDP(5);
t4943 = 4 * MDP(12);
t4942 = 0.2e1 * t4592;
t4941 = pkin(4) / 0.2e1;
t4937 = pkin(4) * t4575;
t4936 = pkin(4) * t4578;
t4935 = pkin(4) * t4581;
t4934 = -qJ(3,1) + qJ(1,1);
t4933 = qJ(3,1) + qJ(1,1);
t4932 = -qJ(3,2) + qJ(1,2);
t4931 = qJ(3,2) + qJ(1,2);
t4930 = -qJ(3,3) + qJ(1,3);
t4929 = qJ(3,3) + qJ(1,3);
t4603 = 0.1e1 / pkin(1);
t4928 = MDP(6) * t4603;
t4927 = MDP(7) * t4603;
t4926 = MDP(8) / pkin(1) ^ 2;
t4925 = qJ(1,1) + 0.2e1 * qJ(2,1);
t4924 = qJ(1,1) - 0.2e1 * qJ(2,1);
t4923 = qJ(1,2) + 0.2e1 * qJ(2,2);
t4922 = qJ(1,2) - 0.2e1 * qJ(2,2);
t4921 = 0.2e1 * qJ(2,3) + qJ(1,3);
t4920 = -0.2e1 * qJ(2,3) + qJ(1,3);
t4919 = MDP(15) * t4603;
t4593 = 0.2e1 * qJ(3,3);
t4452 = (cos(qJ(2,3) - t4930) + cos(qJ(2,3) + t4929)) * t4942 + (-sin(0.2e1 * qJ(3,3) - t4920) + sin(t4593 + t4921) + 0.2e1 * t4576) * pkin(2) + (-sin(qJ(3,3) - t4920) + sin(qJ(3,3) + t4921) + sin(t4930) + sin(t4929)) * pkin(1);
t4568 = qJ(2,3) + qJ(3,3);
t4478 = (-sin(t4593 + qJ(2,3)) + t4575) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - sin(t4568)) * pkin(1);
t4918 = t4452 / t4478;
t4596 = 0.2e1 * qJ(3,2);
t4453 = (cos(qJ(2,2) - t4932) + cos(qJ(2,2) + t4931)) * t4942 + (-sin(0.2e1 * qJ(3,2) - t4922) + sin(t4596 + t4923) + 0.2e1 * t4579) * pkin(2) + (-sin(qJ(3,2) - t4922) + sin(qJ(3,2) + t4923) + sin(t4932) + sin(t4931)) * pkin(1);
t4569 = qJ(2,2) + qJ(3,2);
t4479 = (-sin(t4596 + qJ(2,2)) + t4578) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - sin(t4569)) * pkin(1);
t4917 = t4453 / t4479;
t4599 = 0.2e1 * qJ(3,1);
t4454 = (cos(qJ(2,1) - t4934) + cos(qJ(2,1) + t4933)) * t4942 + (-sin(0.2e1 * qJ(3,1) - t4924) + sin(t4599 + t4925) + 0.2e1 * t4582) * pkin(2) + (-sin(qJ(3,1) - t4924) + sin(qJ(3,1) + t4925) + sin(t4934) + sin(t4933)) * pkin(1);
t4570 = qJ(2,1) + qJ(3,1);
t4480 = (-sin(t4599 + qJ(2,1)) + t4581) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - sin(t4570)) * pkin(1);
t4916 = t4454 / t4480;
t4869 = t4574 * t4584;
t4496 = pkin(2) * t4869 + t4527 * t4575;
t4457 = -t4544 * t4496 - t4953 * t4541;
t4547 = 0.1e1 / t4574;
t4915 = t4457 * t4547;
t4458 = -t4541 * t4496 + t4953 * t4544;
t4914 = t4458 * t4547;
t4865 = t4577 * t4587;
t4497 = pkin(2) * t4865 + t4529 * t4578;
t4459 = -t4545 * t4497 - t4952 * t4542;
t4551 = 0.1e1 / t4577;
t4913 = t4459 * t4551;
t4460 = -t4542 * t4497 + t4952 * t4545;
t4912 = t4460 * t4551;
t4861 = t4580 * t4590;
t4498 = pkin(2) * t4861 + t4531 * t4581;
t4461 = -t4546 * t4498 - t4951 * t4543;
t4555 = 0.1e1 / t4580;
t4911 = t4461 * t4555;
t4462 = -t4543 * t4498 + t4951 * t4546;
t4910 = t4462 * t4555;
t4466 = -t4576 * t4655 + t4585 * t4592;
t4909 = t4466 * t4547;
t4467 = -t4579 * t4654 + t4588 * t4592;
t4908 = t4467 * t4551;
t4468 = -t4582 * t4653 + t4591 * t4592;
t4907 = t4468 * t4555;
t4490 = 0.1e1 / t4655;
t4906 = t4490 * t4547;
t4905 = 0.1e1 / t4655 ^ 2 / t4574 ^ 2;
t4492 = 0.1e1 / t4654;
t4904 = t4492 * t4551;
t4903 = 0.1e1 / t4654 ^ 2 / t4577 ^ 2;
t4494 = 0.1e1 / t4653;
t4902 = t4494 * t4555;
t4901 = 0.1e1 / t4653 ^ 2 / t4580 ^ 2;
t4520 = pkin(2) * cos(t4568) + t4584 * pkin(1);
t4514 = 0.1e1 / t4520;
t4900 = t4514 * t4576;
t4899 = t4514 * t4585;
t4515 = 0.1e1 / t4520 ^ 2;
t4550 = t4576 ^ 2;
t4898 = t4515 * t4550;
t4561 = t4585 ^ 2;
t4897 = t4515 * t4561;
t4521 = pkin(2) * cos(t4569) + t4587 * pkin(1);
t4516 = 0.1e1 / t4521;
t4896 = t4516 * t4579;
t4895 = t4516 * t4588;
t4517 = 0.1e1 / t4521 ^ 2;
t4554 = t4579 ^ 2;
t4894 = t4517 * t4554;
t4564 = t4588 ^ 2;
t4893 = t4517 * t4564;
t4522 = pkin(2) * cos(t4570) + t4590 * pkin(1);
t4518 = 0.1e1 / t4522;
t4892 = t4518 * t4582;
t4891 = t4518 * t4591;
t4519 = 0.1e1 / t4522 ^ 2;
t4558 = t4582 ^ 2;
t4890 = t4519 * t4558;
t4567 = t4591 ^ 2;
t4889 = t4519 * t4567;
t4559 = t4583 ^ 2;
t4523 = pkin(1) * t4583 + t4559 * t4945 - pkin(2);
t4888 = t4523 * t4575;
t4887 = t4523 * t4585;
t4562 = t4586 ^ 2;
t4524 = pkin(1) * t4586 + t4562 * t4945 - pkin(2);
t4886 = t4524 * t4578;
t4885 = t4524 * t4588;
t4565 = t4589 ^ 2;
t4525 = pkin(1) * t4589 + t4565 * t4945 - pkin(2);
t4884 = t4525 * t4581;
t4883 = t4525 * t4591;
t4882 = (pkin(1) + 0.2e1 * t4940) * t4574;
t4881 = (pkin(1) + 0.2e1 * t4939) * t4577;
t4880 = (pkin(1) + 0.2e1 * t4938) * t4580;
t4879 = t4541 * t4576;
t4878 = t4542 * t4579;
t4877 = t4543 * t4582;
t4876 = t4544 * t4576;
t4875 = t4545 * t4579;
t4874 = t4546 * t4582;
t4873 = t4547 * t4585;
t4872 = t4551 * t4588;
t4871 = t4555 * t4591;
t4868 = t4575 * t4584;
t4867 = t4575 * t4585;
t4864 = t4578 * t4587;
t4863 = t4578 * t4588;
t4860 = t4581 * t4590;
t4859 = t4581 * t4591;
t4858 = t4583 * t4574;
t4857 = t4583 * t4584;
t4856 = t4586 * t4577;
t4855 = t4586 * t4587;
t4854 = t4589 * t4580;
t4853 = t4589 * t4590;
t4602 = 0.1e1 / pkin(2);
t4852 = t4602 * t4603;
t4851 = 0.2e1 * t4928;
t4850 = 0.2e1 * t4927;
t4849 = -0.2e1 * t4869;
t4848 = -0.2e1 * t4865;
t4847 = -0.2e1 * t4861;
t4846 = -0.2e1 * t4857;
t4845 = -0.2e1 * t4855;
t4844 = -0.2e1 * t4853;
t4843 = pkin(2) * t4858;
t4842 = pkin(2) * t4856;
t4841 = pkin(2) * t4854;
t4840 = pkin(4) * t4514 * t4584;
t4839 = pkin(4) * t4516 * t4587;
t4838 = pkin(4) * t4518 * t4590;
t4481 = -t4532 * t4870 + (t4559 - 0.1e1) * t4585 * pkin(2);
t4560 = t4584 ^ 2;
t4780 = t4574 * t4867;
t4628 = pkin(1) * t4780 + (t4780 * t4945 + t4532) * t4583;
t4427 = (-t4541 * t4887 + t4544 * t4882) * t4560 + (t4628 * t4541 + t4544 * t4888) * t4584 + t4481 * t4541 - t4544 * t4843;
t4837 = t4427 * t4906;
t4428 = (t4541 * t4882 + t4544 * t4887) * t4560 + (t4541 * t4888 - t4628 * t4544) * t4584 - t4481 * t4544 - t4541 * t4843;
t4836 = t4428 * t4906;
t4482 = -t4533 * t4866 + (t4562 - 0.1e1) * t4588 * pkin(2);
t4563 = t4587 ^ 2;
t4779 = t4577 * t4863;
t4627 = pkin(1) * t4779 + (t4779 * t4945 + t4533) * t4586;
t4429 = (-t4542 * t4885 + t4545 * t4881) * t4563 + (t4627 * t4542 + t4545 * t4886) * t4587 + t4482 * t4542 - t4545 * t4842;
t4835 = t4429 * t4904;
t4430 = (t4542 * t4881 + t4545 * t4885) * t4563 + (t4542 * t4886 - t4627 * t4545) * t4587 - t4482 * t4545 - t4542 * t4842;
t4834 = t4430 * t4904;
t4483 = -t4534 * t4862 + (t4565 - 0.1e1) * t4591 * pkin(2);
t4566 = t4590 ^ 2;
t4778 = t4580 * t4859;
t4626 = pkin(1) * t4778 + (t4778 * t4945 + t4534) * t4589;
t4431 = (-t4543 * t4883 + t4546 * t4880) * t4566 + (t4626 * t4543 + t4546 * t4884) * t4590 + t4483 * t4543 - t4546 * t4841;
t4833 = t4431 * t4902;
t4432 = (t4543 * t4880 + t4546 * t4883) * t4566 + (t4543 * t4884 - t4626 * t4546) * t4590 - t4483 * t4546 - t4543 * t4841;
t4832 = t4432 * t4902;
t4776 = t4918 / 0.2e1;
t4442 = t4603 * t4776;
t4783 = t4547 * t4852;
t4433 = t4466 * t4783 + t4442;
t4831 = t4433 * t4906;
t4775 = t4917 / 0.2e1;
t4443 = t4603 * t4775;
t4782 = t4551 * t4852;
t4434 = t4467 * t4782 + t4443;
t4830 = t4434 * t4904;
t4774 = t4916 / 0.2e1;
t4444 = t4603 * t4774;
t4781 = t4555 * t4852;
t4435 = t4468 * t4781 + t4444;
t4829 = t4435 * t4902;
t4828 = t4490 * t4873;
t4827 = t4603 * t4906;
t4826 = t4492 * t4872;
t4825 = t4603 * t4904;
t4824 = t4494 * t4871;
t4823 = t4603 * t4902;
t4508 = t4857 - t4870;
t4822 = t4508 * t4899;
t4509 = t4855 - t4866;
t4821 = t4509 * t4895;
t4510 = t4853 - t4862;
t4820 = t4510 * t4891;
t4511 = t4575 * t4583 + t4869;
t4819 = t4511 * t4899;
t4512 = t4578 * t4586 + t4865;
t4818 = t4512 * t4895;
t4513 = t4581 * t4589 + t4861;
t4817 = t4513 * t4891;
t4816 = t4514 * t4879;
t4815 = t4514 * t4876;
t4814 = t4547 * t4900;
t4813 = t4514 * t4873;
t4812 = t4575 * t4900;
t4811 = t4514 * t4867;
t4535 = t4541 ^ 2;
t4810 = t4535 * t4898;
t4538 = t4544 ^ 2;
t4809 = t4538 * t4898;
t4549 = t4575 ^ 2;
t4808 = t4549 * t4898;
t4807 = t4515 * t4868;
t4806 = t4515 * t4576 * t4585;
t4805 = t4516 * t4878;
t4804 = t4516 * t4875;
t4803 = t4551 * t4896;
t4802 = t4516 * t4872;
t4801 = t4578 * t4896;
t4800 = t4516 * t4863;
t4536 = t4542 ^ 2;
t4799 = t4536 * t4894;
t4539 = t4545 ^ 2;
t4798 = t4539 * t4894;
t4553 = t4578 ^ 2;
t4797 = t4553 * t4894;
t4796 = t4517 * t4864;
t4795 = t4517 * t4579 * t4588;
t4794 = t4518 * t4877;
t4793 = t4518 * t4874;
t4792 = t4555 * t4892;
t4791 = t4518 * t4871;
t4790 = t4581 * t4892;
t4789 = t4518 * t4859;
t4537 = t4543 ^ 2;
t4788 = t4537 * t4890;
t4540 = t4546 ^ 2;
t4787 = t4540 * t4890;
t4557 = t4581 ^ 2;
t4786 = t4557 * t4890;
t4785 = t4519 * t4860;
t4784 = t4519 * t4582 * t4591;
t4410 = t4430 * t4825;
t4379 = t4460 * t4782 - t4410;
t4777 = t4379 * t4895;
t4773 = t4852 / 0.2e1;
t4772 = t4514 * t4560 * t4946;
t4771 = t4516 * t4563 * t4946;
t4770 = t4518 * t4566 * t4946;
t4769 = t4576 * t4840;
t4768 = t4585 * t4840;
t4767 = t4579 * t4839;
t4766 = t4588 * t4839;
t4765 = t4582 * t4838;
t4764 = t4591 * t4838;
t4763 = pkin(4) * t4811;
t4762 = pkin(4) * t4800;
t4761 = pkin(4) * t4789;
t4760 = t4899 * t4918;
t4759 = t4895 * t4917;
t4758 = t4891 * t4916;
t4757 = t4508 * t4816;
t4756 = t4508 * t4815;
t4755 = t4508 * t4813;
t4754 = t4509 * t4805;
t4753 = t4509 * t4804;
t4752 = t4509 * t4802;
t4751 = t4510 * t4794;
t4750 = t4510 * t4793;
t4749 = t4510 * t4791;
t4748 = t4511 * t4816;
t4747 = t4511 * t4815;
t4746 = t4511 * t4813;
t4745 = t4512 * t4805;
t4744 = t4512 * t4804;
t4743 = t4512 * t4802;
t4742 = t4513 * t4794;
t4741 = t4513 * t4793;
t4740 = t4513 * t4791;
t4739 = t4541 * t4814;
t4738 = t4541 * t4812;
t4737 = t4544 * t4814;
t4736 = t4544 * t4812;
t4735 = t4898 * t4950;
t4734 = t4541 * t4806;
t4733 = t4544 * t4806;
t4732 = t4549 * t4806;
t4731 = t4550 * t4807;
t4730 = t4542 * t4803;
t4729 = t4542 * t4801;
t4728 = t4545 * t4803;
t4727 = t4545 * t4801;
t4726 = t4894 * t4949;
t4725 = t4542 * t4795;
t4724 = t4545 * t4795;
t4723 = t4553 * t4795;
t4722 = t4554 * t4796;
t4721 = t4543 * t4792;
t4720 = t4543 * t4790;
t4719 = t4546 * t4792;
t4718 = t4546 * t4790;
t4717 = t4890 * t4948;
t4716 = t4543 * t4784;
t4715 = t4546 * t4784;
t4714 = t4557 * t4784;
t4713 = t4558 * t4785;
t4712 = t4547 * t4773;
t4711 = t4551 * t4773;
t4710 = t4555 * t4773;
t4709 = t4576 * t4772;
t4708 = t4579 * t4771;
t4707 = t4582 * t4770;
t4706 = t4574 * t4769;
t4705 = t4583 * t4769;
t4704 = t4577 * t4767;
t4703 = t4586 * t4767;
t4702 = t4580 * t4765;
t4701 = t4589 * t4765;
t4700 = pkin(4) * t4738;
t4699 = pkin(4) * t4736;
t4698 = pkin(4) * t4729;
t4697 = pkin(4) * t4727;
t4696 = pkin(4) * t4720;
t4695 = pkin(4) * t4718;
t4694 = t4490 * t4755;
t4693 = t4490 * t4746;
t4692 = t4492 * t4752;
t4691 = t4492 * t4743;
t4690 = t4494 * t4749;
t4689 = t4494 * t4740;
t4688 = t4508 * t4739;
t4687 = t4508 * t4737;
t4686 = t4509 * t4730;
t4685 = t4509 * t4728;
t4684 = t4510 * t4721;
t4683 = t4510 * t4719;
t4682 = t4511 * t4739;
t4681 = t4511 * t4737;
t4680 = t4512 * t4730;
t4679 = t4512 * t4728;
t4678 = t4513 * t4721;
t4677 = t4513 * t4719;
t4676 = t4806 * t4868;
t4675 = t4795 * t4864;
t4674 = t4784 * t4860;
t4673 = t4776 * t4906;
t4672 = -t4760 / 0.2e1;
t4671 = t4775 * t4904;
t4670 = -t4759 / 0.2e1;
t4669 = t4774 * t4902;
t4668 = -t4758 / 0.2e1;
t4667 = t4776 * t4879;
t4666 = t4775 * t4878;
t4665 = t4774 * t4877;
t4664 = -t4876 * t4918 / 0.2e1;
t4663 = -t4875 * t4917 / 0.2e1;
t4662 = -t4874 * t4916 / 0.2e1;
t4661 = t4541 * t4706;
t4660 = t4541 * t4705;
t4659 = t4542 * t4704;
t4658 = t4542 * t4703;
t4657 = t4543 * t4702;
t4656 = t4543 * t4701;
t4652 = t4427 * t4490 * t4739;
t4651 = t4428 * t4490 * t4737;
t4650 = t4429 * t4492 * t4730;
t4649 = t4430 * t4492 * t4728;
t4648 = t4431 * t4494 * t4721;
t4647 = t4432 * t4494 * t4719;
t4646 = t4490 * t4688;
t4645 = t4490 * t4687;
t4644 = t4490 * t4682;
t4643 = t4490 * t4681;
t4642 = t4492 * t4686;
t4641 = t4492 * t4685;
t4640 = t4492 * t4680;
t4639 = t4492 * t4679;
t4638 = t4494 * t4684;
t4637 = t4494 * t4683;
t4636 = t4494 * t4678;
t4635 = t4494 * t4677;
t4634 = t4514 * t4667;
t4633 = t4514 * t4664;
t4632 = t4516 * t4666;
t4631 = t4516 * t4663;
t4630 = t4518 * t4665;
t4629 = t4518 * t4662;
t4625 = t4433 * t4937 + t4585 * t4772;
t4624 = t4434 * t4936 + t4588 * t4771;
t4623 = t4435 * t4935 + t4591 * t4770;
t4463 = (t4559 - 0.1e1 / 0.2e1) * t4868 + (t4560 - 0.1e1 / 0.2e1) * t4858;
t4464 = (t4562 - 0.1e1 / 0.2e1) * t4864 + (t4563 - 0.1e1 / 0.2e1) * t4856;
t4465 = (t4565 - 0.1e1 / 0.2e1) * t4860 + (t4566 - 0.1e1 / 0.2e1) * t4854;
t4505 = t4511 ^ 2;
t4506 = t4512 ^ 2;
t4507 = t4513 ^ 2;
t4611 = t4494 * (-t4431 * t4546 + t4432 * t4543) * t4792;
t4612 = t4492 * (-t4429 * t4545 + t4430 * t4542) * t4803;
t4613 = t4490 * (-t4427 * t4544 + t4428 * t4541) * t4814;
t4622 = (-t4575 * t4613 - t4578 * t4612 - t4581 * t4611) * t4928 + (-t4584 * t4613 - t4587 * t4612 - t4590 * t4611) * t4927 + (t4427 * t4428 * t4905 + t4429 * t4430 * t4903 + t4431 * t4432 * t4901) * t4926 + (-t4463 * t4735 - t4464 * t4726 - t4465 * t4717) * t4943 + (-t4505 * t4735 - t4506 * t4726 - t4507 * t4717) * MDP(11) + (-t4713 * t4948 - t4722 * t4949 - t4731 * t4950) * t4944 + (-t4549 * t4735 - t4553 * t4726 - t4557 * t4717) * MDP(4) + (-t4717 - t4726 - t4735) * MDP(1);
t4606 = t4518 * (t4431 * t4824 + t4665);
t4608 = t4516 * (t4429 * t4826 + t4666);
t4610 = t4514 * (t4427 * t4828 + t4667);
t4621 = (t4575 * t4610 + t4578 * t4608 + t4581 * t4606) * t4928 + (t4584 * t4610 + t4587 * t4608 + t4590 * t4606) * t4927 + (-t4427 * t4673 - t4429 * t4671 - t4431 * t4669) * t4926 + (-t4463 * t4734 - t4464 * t4725 - t4465 * t4716) * t4943 + (-t4505 * t4734 - t4506 * t4725 - t4507 * t4716) * MDP(11) + (-t4541 * t4676 - t4542 * t4675 - t4543 * t4674) * t4944 + (-t4541 * t4732 - t4542 * t4723 - t4543 * t4714) * MDP(4) + (-t4716 - t4725 - t4734) * MDP(1);
t4605 = t4518 * (t4432 * t4824 + t4662);
t4607 = t4516 * (t4430 * t4826 + t4663);
t4609 = t4514 * (t4428 * t4828 + t4664);
t4620 = (t4575 * t4609 + t4578 * t4607 + t4581 * t4605) * t4928 + (t4584 * t4609 + t4587 * t4607 + t4590 * t4605) * t4927 + (-t4428 * t4673 - t4430 * t4671 - t4432 * t4669) * t4926 + (t4463 * t4733 + t4464 * t4724 + t4465 * t4715) * t4943 + (t4505 * t4733 + t4506 * t4724 + t4507 * t4715) * MDP(11) + (t4544 * t4676 + t4545 * t4675 + t4546 * t4674) * t4944 + (t4544 * t4732 + t4545 * t4723 + t4546 * t4714) * MDP(4) + (t4715 + t4724 + t4733) * MDP(1);
t4407 = t4427 * t4827;
t4376 = t4457 * t4783 - t4407;
t4619 = -t4376 * t4937 + t4541 * t4709;
t4408 = t4428 * t4827;
t4377 = t4458 * t4783 - t4408;
t4618 = t4377 * t4937 + t4544 * t4709;
t4409 = t4429 * t4825;
t4378 = t4459 * t4782 - t4409;
t4617 = -t4378 * t4936 + t4542 * t4708;
t4616 = t4379 * t4936 + t4545 * t4708;
t4411 = t4431 * t4823;
t4380 = t4461 * t4781 - t4411;
t4615 = -t4380 * t4935 + t4543 * t4707;
t4412 = t4432 * t4823;
t4381 = t4462 * t4781 - t4412;
t4614 = t4381 * t4935 + t4546 * t4707;
t4489 = t4589 * t4764;
t4488 = t4586 * t4766;
t4487 = t4583 * t4768;
t4486 = t4580 * t4764;
t4485 = t4577 * t4766;
t4484 = t4574 * t4768;
t4474 = t4546 * t4701;
t4473 = t4545 * t4703;
t4472 = t4544 * t4705;
t4471 = t4546 * t4702;
t4470 = t4545 * t4704;
t4469 = t4544 * t4706;
t4438 = t4761 + t4774;
t4437 = t4762 + t4775;
t4436 = t4763 + t4776;
t4426 = -t4438 * t4580 + t4489;
t4425 = -t4437 * t4577 + t4488;
t4424 = -t4436 * t4574 + t4487;
t4423 = t4438 * t4589 + t4486;
t4422 = t4437 * t4586 + t4485;
t4421 = t4436 * t4583 + t4484;
t4420 = -pkin(1) * t4789 + t4435 * t4941;
t4419 = t4761 + (t4468 * t4710 + t4444) * t4946;
t4418 = -pkin(1) * t4800 + t4434 * t4941;
t4417 = t4762 + (t4467 * t4711 + t4443) * t4946;
t4416 = -pkin(1) * t4811 + t4433 * t4941;
t4415 = t4763 + (t4466 * t4712 + t4442) * t4946;
t4405 = -t4419 * t4580 + t4489;
t4404 = t4419 * t4589 + t4486;
t4403 = -t4417 * t4577 + t4488;
t4402 = t4417 * t4586 + t4485;
t4401 = -t4415 * t4574 + t4487;
t4400 = t4415 * t4583 + t4484;
t4399 = t4695 - t4832;
t4398 = -t4696 - t4833;
t4397 = t4697 - t4834;
t4396 = -t4698 - t4835;
t4395 = t4699 - t4836;
t4394 = -t4700 - t4837;
t4393 = -t4399 * t4580 + t4474;
t4392 = -t4398 * t4580 - t4656;
t4391 = -t4397 * t4577 + t4473;
t4390 = -t4396 * t4577 - t4658;
t4389 = -t4395 * t4574 + t4472;
t4388 = -t4394 * t4574 - t4660;
t4387 = t4399 * t4589 + t4471;
t4386 = t4398 * t4589 - t4657;
t4385 = t4397 * t4586 + t4470;
t4384 = t4396 * t4586 - t4659;
t4383 = t4395 * t4583 + t4469;
t4382 = t4394 * t4583 - t4661;
t4375 = -pkin(1) * t4718 + t4381 * t4941;
t4374 = pkin(1) * t4720 + t4380 * t4941;
t4373 = t4695 + (t4462 * t4710 - t4412) * t4946;
t4372 = t4696 + (t4461 * t4710 - t4411) * t4947;
t4371 = -pkin(1) * t4727 + t4379 * t4941;
t4370 = pkin(1) * t4729 + t4378 * t4941;
t4369 = t4697 + (t4460 * t4711 - t4410) * t4946;
t4368 = t4698 + (t4459 * t4711 - t4409) * t4947;
t4367 = -pkin(1) * t4736 + t4377 * t4941;
t4366 = pkin(1) * t4738 + t4376 * t4941;
t4365 = t4699 + (t4458 * t4712 - t4408) * t4946;
t4364 = t4700 + (t4457 * t4712 - t4407) * t4947;
t4363 = -t4373 * t4580 + t4474;
t4362 = t4373 * t4589 + t4471;
t4361 = -t4372 * t4589 - t4657;
t4360 = t4372 * t4580 - t4656;
t4359 = -t4369 * t4577 + t4473;
t4358 = t4369 * t4586 + t4470;
t4357 = -t4368 * t4586 - t4659;
t4356 = t4368 * t4577 - t4658;
t4355 = -t4365 * t4574 + t4472;
t4354 = t4365 * t4583 + t4469;
t4353 = -t4364 * t4583 - t4661;
t4352 = t4364 * t4574 - t4660;
t4351 = t4420 * t4844 + t4623 * t4580;
t4350 = t4420 * t4847 - t4623 * t4589;
t4349 = t4418 * t4845 + t4624 * t4577;
t4348 = t4418 * t4848 - t4624 * t4586;
t4347 = t4416 * t4846 + t4625 * t4574;
t4346 = t4416 * t4849 - t4625 * t4583;
t4345 = t4375 * t4844 + t4614 * t4580;
t4344 = t4374 * t4844 - t4615 * t4580;
t4343 = t4375 * t4847 - t4614 * t4589;
t4342 = t4374 * t4847 + t4615 * t4589;
t4341 = t4371 * t4845 + t4616 * t4577;
t4340 = t4370 * t4845 - t4617 * t4577;
t4339 = t4371 * t4848 - t4616 * t4586;
t4338 = t4370 * t4848 + t4617 * t4586;
t4337 = t4367 * t4846 + t4618 * t4574;
t4336 = t4366 * t4846 - t4619 * t4574;
t4335 = t4367 * t4849 - t4618 * t4583;
t4334 = t4366 * t4849 + t4619 * t4583;
t1 = [(t4787 + t4798 + t4809) * MDP(1) + (t4538 * t4808 + t4539 * t4797 + t4540 * t4786) * MDP(4) + (t4538 * t4731 + t4539 * t4722 + t4540 * t4713) * t4944 + (t4575 * t4651 + t4578 * t4649 + t4581 * t4647) * t4851 + (t4584 * t4651 + t4587 * t4649 + t4590 * t4647) * t4850 + (t4428 ^ 2 * t4905 + t4430 ^ 2 * t4903 + t4432 ^ 2 * t4901) * t4926 + (t4505 * t4809 + t4506 * t4798 + t4507 * t4787) * MDP(11) + (t4463 * t4809 + t4464 * t4798 + t4465 * t4787) * t4943 + (-t4377 * t4747 - t4379 * t4744 - t4381 * t4741 + (t4428 * t4643 + t4430 * t4639 + t4432 * t4635 + (-t4458 * t4681 - t4460 * t4679 - t4462 * t4677) * t4602) * t4603) * MDP(13) + (-t4377 * t4756 - t4379 * t4753 - t4381 * t4750 + (t4428 * t4645 + t4430 * t4641 + t4432 * t4637 + (-t4458 * t4687 - t4460 * t4685 - t4462 * t4683) * t4602) * t4603) * MDP(14) + (-t4377 * t4836 - t4379 * t4834 - t4381 * t4832 + (t4377 * t4914 + t4379 * t4912 + t4381 * t4910) * t4602) * t4919 + (-t4335 * t4815 - t4339 * t4804 - t4343 * t4793 + (-t4354 * t4836 - t4358 * t4834 - t4362 * t4832 + (t4383 * t4914 + t4385 * t4912 + t4387 * t4910) * t4602) * t4603) * MDP(16) + (-t4337 * t4815 - t4341 * t4804 - t4345 * t4793 + (-t4355 * t4836 - t4359 * t4834 - t4363 * t4832 + (t4389 * t4914 + t4391 * t4912 + t4393 * t4910) * t4602) * t4603) * MDP(17) + MDP(18); (-t4376 * t4747 - t4378 * t4744 - t4380 * t4741 + (-t4428 * t4644 - t4430 * t4640 - t4432 * t4636 + (t4458 * t4682 + t4460 * t4680 + t4462 * t4678) * t4602) * t4603) * MDP(13) + (-t4376 * t4756 - t4378 * t4753 - t4380 * t4750 + (-t4428 * t4646 - t4430 * t4642 - t4432 * t4638 + (t4458 * t4688 + t4460 * t4686 + t4462 * t4684) * t4602) * t4603) * MDP(14) + (-t4376 * t4836 - t4378 * t4834 - t4380 * t4832 + (t4376 * t4914 + t4378 * t4912 + t4380 * t4910) * t4602) * t4919 + (-t4334 * t4815 - t4338 * t4804 - t4342 * t4793 + (-t4353 * t4836 - t4357 * t4834 - t4361 * t4832 + (t4382 * t4914 + t4384 * t4912 + t4386 * t4910) * t4602) * t4603) * MDP(16) + (-t4336 * t4815 - t4340 * t4804 - t4344 * t4793 + (-t4352 * t4836 - t4356 * t4834 - t4360 * t4832 + (t4388 * t4914 + t4390 * t4912 + t4392 * t4910) * t4602) * t4603) * MDP(17) + t4622; (-t4433 * t4747 - t4434 * t4744 - t4435 * t4741 + (t4428 * t4693 + t4430 * t4691 + t4432 * t4689 + (-t4458 * t4746 - t4460 * t4743 - t4462 * t4740) * t4602) * t4603) * MDP(13) + (-t4433 * t4756 - t4434 * t4753 - t4435 * t4750 + (t4428 * t4694 + t4430 * t4692 + t4432 * t4690 + (-t4458 * t4755 - t4460 * t4752 - t4462 * t4749) * t4602) * t4603) * MDP(14) + (-t4428 * t4831 - t4430 * t4830 - t4432 * t4829 + (t4433 * t4914 + t4434 * t4912 + t4435 * t4910) * t4602) * t4919 + (-t4346 * t4815 - t4348 * t4804 - t4350 * t4793 + (-t4400 * t4836 - t4402 * t4834 - t4404 * t4832 + (t4421 * t4914 + t4422 * t4912 + t4423 * t4910) * t4602) * t4603) * MDP(16) + (-t4347 * t4815 - t4349 * t4804 - t4351 * t4793 + (-t4401 * t4836 - t4403 * t4834 - t4405 * t4832 + (t4424 * t4914 + t4425 * t4912 + t4426 * t4910) * t4602) * t4603) * MDP(17) + t4620; (t4377 * t4748 + t4379 * t4745 + t4381 * t4742 + (t4427 * t4643 + t4429 * t4639 + t4431 * t4635 + (-t4457 * t4681 - t4459 * t4679 - t4461 * t4677) * t4602) * t4603) * MDP(13) + (t4377 * t4757 + t4379 * t4754 + t4381 * t4751 + (t4427 * t4645 + t4429 * t4641 + t4431 * t4637 + (-t4457 * t4687 - t4459 * t4685 - t4461 * t4683) * t4602) * t4603) * MDP(14) + (-t4377 * t4837 - t4379 * t4835 - t4381 * t4833 + (t4377 * t4915 + t4379 * t4913 + t4381 * t4911) * t4602) * t4919 + (t4335 * t4816 + t4339 * t4805 + t4343 * t4794 + (-t4354 * t4837 - t4358 * t4835 - t4362 * t4833 + (t4383 * t4915 + t4385 * t4913 + t4387 * t4911) * t4602) * t4603) * MDP(16) + (t4337 * t4816 + t4341 * t4805 + t4345 * t4794 + (-t4355 * t4837 - t4359 * t4835 - t4363 * t4833 + (t4389 * t4915 + t4391 * t4913 + t4393 * t4911) * t4602) * t4603) * MDP(17) + t4622; (t4788 + t4799 + t4810) * MDP(1) + (t4535 * t4808 + t4536 * t4797 + t4537 * t4786) * MDP(4) + (t4535 * t4731 + t4536 * t4722 + t4537 * t4713) * t4944 + (-t4575 * t4652 - t4578 * t4650 - t4581 * t4648) * t4851 + (-t4584 * t4652 - t4587 * t4650 - t4590 * t4648) * t4850 + (t4427 ^ 2 * t4905 + t4429 ^ 2 * t4903 + t4431 ^ 2 * t4901) * t4926 + (t4505 * t4810 + t4506 * t4799 + t4507 * t4788) * MDP(11) + (t4463 * t4810 + t4464 * t4799 + t4465 * t4788) * t4943 + (t4376 * t4748 + t4378 * t4745 + t4380 * t4742 + (-t4427 * t4644 - t4429 * t4640 - t4431 * t4636 + (t4457 * t4682 + t4459 * t4680 + t4461 * t4678) * t4602) * t4603) * MDP(13) + (t4376 * t4757 + t4378 * t4754 + t4380 * t4751 + (-t4427 * t4646 - t4429 * t4642 - t4431 * t4638 + (t4457 * t4688 + t4459 * t4686 + t4461 * t4684) * t4602) * t4603) * MDP(14) + (-t4376 * t4837 - t4378 * t4835 - t4380 * t4833 + (t4376 * t4915 + t4378 * t4913 + t4380 * t4911) * t4602) * t4919 + (t4334 * t4816 + t4338 * t4805 + t4342 * t4794 + (-t4353 * t4837 - t4357 * t4835 - t4361 * t4833 + (t4382 * t4915 + t4384 * t4913 + t4386 * t4911) * t4602) * t4603) * MDP(16) + (t4336 * t4816 + t4340 * t4805 + t4344 * t4794 + (-t4352 * t4837 - t4356 * t4835 - t4360 * t4833 + (t4388 * t4915 + t4390 * t4913 + t4392 * t4911) * t4602) * t4603) * MDP(17) + MDP(18); (t4433 * t4748 + t4434 * t4745 + t4435 * t4742 + (t4427 * t4693 + t4429 * t4691 + t4431 * t4689 + (-t4457 * t4746 - t4459 * t4743 - t4461 * t4740) * t4602) * t4603) * MDP(13) + (t4433 * t4757 + t4434 * t4754 + t4435 * t4751 + (t4427 * t4694 + t4429 * t4692 + t4431 * t4690 + (-t4457 * t4755 - t4459 * t4752 - t4461 * t4749) * t4602) * t4603) * MDP(14) + (-t4427 * t4831 - t4429 * t4830 - t4431 * t4829 + (t4433 * t4915 + t4434 * t4913 + t4435 * t4911) * t4602) * t4919 + (t4346 * t4816 + t4348 * t4805 + t4350 * t4794 + (-t4400 * t4837 - t4402 * t4835 - t4404 * t4833 + (t4421 * t4915 + t4422 * t4913 + t4423 * t4911) * t4602) * t4603) * MDP(16) + (t4347 * t4816 + t4349 * t4805 + t4351 * t4794 + (-t4401 * t4837 - t4403 * t4835 - t4405 * t4833 + (t4424 * t4915 + t4425 * t4913 + t4426 * t4911) * t4602) * t4603) * MDP(17) + t4621; (-t4377 * t4819 - t4512 * t4777 - t4381 * t4817 + (t4513 * t4629 + t4512 * t4631 + t4511 * t4633 + (-t4466 * t4681 - t4467 * t4679 - t4468 * t4677) * t4602) * t4603) * MDP(13) + (-t4377 * t4822 - t4509 * t4777 - t4381 * t4820 + (t4510 * t4629 + t4509 * t4631 + t4508 * t4633 + (-t4466 * t4687 - t4467 * t4685 - t4468 * t4683) * t4602) * t4603) * MDP(14) + (t4381 * t4774 + t4379 * t4775 + t4377 * t4776 + (t4377 * t4909 + t4379 * t4908 + t4381 * t4907) * t4602) * t4919 + (-t4335 * t4899 - t4339 * t4895 - t4343 * t4891 + (t4362 * t4774 + t4358 * t4775 + t4354 * t4776 + (t4383 * t4909 + t4385 * t4908 + t4387 * t4907) * t4602) * t4603) * MDP(16) + (-t4337 * t4899 - t4341 * t4895 - t4345 * t4891 + (t4363 * t4774 + t4359 * t4775 + t4355 * t4776 + (t4389 * t4909 + t4391 * t4908 + t4393 * t4907) * t4602) * t4603) * MDP(17) + t4620; (-t4376 * t4819 - t4378 * t4818 - t4380 * t4817 + (t4513 * t4630 + t4512 * t4632 + t4511 * t4634 + (t4466 * t4682 + t4467 * t4680 + t4468 * t4678) * t4602) * t4603) * MDP(13) + (-t4376 * t4822 - t4378 * t4821 - t4380 * t4820 + (t4510 * t4630 + t4509 * t4632 + t4508 * t4634 + (t4466 * t4688 + t4467 * t4686 + t4468 * t4684) * t4602) * t4603) * MDP(14) + (t4380 * t4774 + t4378 * t4775 + t4376 * t4776 + (t4376 * t4909 + t4378 * t4908 + t4380 * t4907) * t4602) * t4919 + (-t4334 * t4899 - t4338 * t4895 - t4342 * t4891 + (t4361 * t4774 + t4357 * t4775 + t4353 * t4776 + (t4382 * t4909 + t4384 * t4908 + t4386 * t4907) * t4602) * t4603) * MDP(16) + (-t4336 * t4899 - t4340 * t4895 - t4344 * t4891 + (t4360 * t4774 + t4356 * t4775 + t4352 * t4776 + (t4388 * t4909 + t4390 * t4908 + t4392 * t4907) * t4602) * t4603) * MDP(17) + t4621; (t4889 + t4893 + t4897) * MDP(1) + (t4549 * t4897 + t4553 * t4893 + t4557 * t4889) * MDP(4) + (t4561 * t4807 + t4564 * t4796 + t4567 * t4785) * t4944 + (-t4575 * t4760 - t4578 * t4759 - t4581 * t4758) * t4928 + (-t4584 * t4760 - t4587 * t4759 - t4590 * t4758) * t4927 + (t4454 ^ 2 / t4480 ^ 2 / 0.4e1 + t4453 ^ 2 / t4479 ^ 2 / 0.4e1 + t4452 ^ 2 / t4478 ^ 2 / 0.4e1) * t4926 + (t4505 * t4897 + t4506 * t4893 + t4507 * t4889) * MDP(11) + (t4463 * t4897 + t4464 * t4893 + t4465 * t4889) * t4943 + (-t4433 * t4819 - t4434 * t4818 - t4435 * t4817 + (t4513 * t4668 + t4512 * t4670 + t4511 * t4672 + (-t4466 * t4746 - t4467 * t4743 - t4468 * t4740) * t4602) * t4603) * MDP(13) + (-t4433 * t4822 - t4434 * t4821 - t4435 * t4820 + (t4510 * t4668 + t4509 * t4670 + t4508 * t4672 + (-t4466 * t4755 - t4467 * t4752 - t4468 * t4749) * t4602) * t4603) * MDP(14) + (t4435 * t4774 + t4434 * t4775 + t4433 * t4776 + (t4433 * t4909 + t4434 * t4908 + t4435 * t4907) * t4602) * t4919 + (-t4346 * t4899 - t4348 * t4895 - t4350 * t4891 + (t4404 * t4774 + t4402 * t4775 + t4400 * t4776 + (t4421 * t4909 + t4422 * t4908 + t4423 * t4907) * t4602) * t4603) * MDP(16) + (-t4347 * t4899 - t4349 * t4895 - t4351 * t4891 + (t4405 * t4774 + t4403 * t4775 + t4401 * t4776 + (t4424 * t4909 + t4425 * t4908 + t4426 * t4907) * t4602) * t4603) * MDP(17) + MDP(18);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
