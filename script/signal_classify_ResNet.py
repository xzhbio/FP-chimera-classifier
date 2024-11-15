#导入相关模块
from torch.utils.data import DataLoader,Dataset
import argparse
import os
import torch 
from tqdm import tqdm
# from torchvision import transforms
import numpy as np
from torch import nn

SEED = 1234

np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)
torch.backends.cudnn.deterministic = True

class Signal_Data(Dataset): #继承Dataset
    def __init__(self, data_path): #__init__是初始化该类的一些基础参数
        self.data = np.load(data_path, allow_pickle=True)
            
            # self.data = [x for x in data if x[1]==1]
    def __len__(self):#返回整个数据集的大小
        return len(self.data)
    
    def __getitem__(self,index):#根据索引index返回dataset[index]
        signal = self.data[index][0].astype(np.float32)
        # signal = (signal - np.mean(signal)) / np.std(signal)
        label = self.data[index][1]
        reads_id = self.data[index][2]

        return signal,label,reads_id#返回该样本


class ResidualBlock(nn.Module):
    expansion = 1 #每一个conv的卷积核个数的倍数
    def __init__(self, in_channel, out_channel, stride=1, downsample=None):#downsample对应虚线残差结构
        super(ResidualBlock, self).__init__()
        self.conv1 = nn.Conv1d(in_channels=in_channel, out_channels=out_channel,
                               kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm1d(out_channel)#BN处理
        self.relu = nn.ReLU()
        self.conv2 = nn.Conv1d(in_channels=out_channel, out_channels=out_channel,
                               kernel_size=3, stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm1d(out_channel)
        self.downsample = downsample

    def forward(self, x):
        identity = x #捷径上的输出值
        if self.downsample is not None:
            identity = self.downsample(x)

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)

        out += identity
        out = self.relu(out)

        return out


# 创建一个包含多个残差块的模型
class ResNet(nn.Module):
    def __init__(self, num_classes=2, include_top=True):
        super(ResNet, self).__init__()
        self.include_top = include_top
        self.in_channel = 64

        self.conv1 = nn.Conv1d(1, 64, kernel_size=7, stride=2, padding=3)
        self.bn1 = nn.BatchNorm1d(64)
        self.relu = nn.ReLU(inplace=True)
        self.maxpool = nn.MaxPool1d(kernel_size=3, stride=2, padding=1)

        self.maxpool_cp = nn.MaxPool1d(kernel_size=7, stride=2, padding=1)
        self.block1 = self._make_layer(ResidualBlock, 64, 2)
        self.block2 = self._make_layer(ResidualBlock, 128, 2, stride=2)
        self.block3 = self._make_layer(ResidualBlock, 256, 2, stride=2)
        self.block4 = self._make_layer(ResidualBlock, 512, 2, stride=2)

        self.avgpool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Linear(512, num_classes)

    def _make_layer(self, block, channel, block_num, stride=1):
        downsample = None
        if stride != 1 or self.in_channel != channel * block.expansion:
            downsample = nn.Sequential(
                nn.Conv1d(self.in_channel, channel * block.expansion, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm1d(channel * block.expansion))

        layers = []
        layers.append(block(self.in_channel, channel, downsample=downsample, stride=stride))
        self.in_channel = channel * block.expansion

        for _ in range(1, block_num):
            layers.append(block(self.in_channel, channel))

        return nn.Sequential(*layers)


    def forward(self, x):
        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.maxpool(out)
        # out = self.cnn(x)

        out = self.maxpool_cp(self.block1(out))
        out = self.maxpool_cp(self.block2(out))
        out = self.maxpool_cp(self.block3(out))
        out = self.maxpool_cp(self.block4(out))

        out = self.avgpool(out)
        out = torch.flatten(out, 1)
        out = self.fc(out)
        return out
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--train',type=str)
    parser.add_argument('--external',type=str, help="External Validation Set")
    parser.add_argument('--save',type=str)
    parser.add_argument('--load',type=str)
    parser.add_argument('--epoch',type=int,default=20)
    # parser.add_argument('--roc',type=str, required=False)
    parser.add_argument('--cuda',type=str,default='0')
    parser.add_argument('--lr',type=float,default=1e-5)
    parser.add_argument('--wd',type=float,default=1e-2)
    parser.add_argument('--batch_size',type=int,default=512)
    parser.add_argument('--output','-o',type=str)
    args = parser.parse_args()

    os.environ['CUDA_VISIBLE_DEVICES'] = args.cuda
    torch.set_num_threads(5)
    if args.train:
        all_data = Signal_Data(args.train)
        n_train = int(len(all_data) * 0.8)
        n_valid  = len(all_data) - n_train
        train_data, temp_data = torch.utils.data.random_split(all_data,[n_train,n_valid])
        test_data, valid_data = torch.utils.data.random_split(temp_data,[int(len(temp_data) * 0.5),len(temp_data) - int(len(temp_data) * 0.5)])

        # test1_data, valid_data = torch.utils.data.random_split(test_data,[int(len(test_data) * 0.5),len(test_data) - int(len(test_data) * 0.5)])
        test_dataloader = DataLoader(test_data, batch_size=args.batch_size, shuffle=True)
        train_dataloader = DataLoader(train_data, batch_size=args.batch_size, shuffle=True)
        # test1_dataloader = DataLoader(test1_data, batch_size=120, shuffle=True)
        valid_dataloader = DataLoader(valid_data, batch_size=args.batch_size, shuffle=True)

    if args.external:
        external_data = Signal_Data(args.external)
        external_dataloader = DataLoader(external_data, batch_size=args.batch_size, shuffle=True)

    model = ResNet(num_classes=2)

    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.wd)
    model.cuda()
    model = nn.DataParallel(model)

    if args.external:
        model.load_state_dict(torch.load(args.load))
        model.eval()
    else:
        model.train()

    num_epochs = args.epoch
    max_acc = 0
    if args.train:
        for epoch in range(num_epochs):
            total = 0
            correct = 0
            total_loss = 0
            total_batches = len(train_dataloader)
            if not args.load:
                with tqdm(total=total_batches) as pbar:
                    for batch_idx, (signals, labels, _) in enumerate(train_dataloader):
                        signals = signals.reshape(signals.shape[0], 1, signals.shape[1])
                        # print(signals.shape)
                        signals = signals.cuda()
                        labels = labels.cuda()
                        outputs = model(signals)
                        _, predicted = torch.max(outputs.data, 1)
                        correct += (predicted == labels).sum().item()
                        total += labels.size(0)
                        loss = criterion(outputs, labels)
                        total_loss += loss.item() * signals.size(0)
                        # print(loss)
                        optimizer.zero_grad()
                        loss.backward()
                        optimizer.step()
                        pbar.update(1)
                    accuracy = correct / total
                    average_loss = total_loss / len(train_data)
                    print(f"Train Accuracy: {accuracy}; Training Loss: {average_loss}; epoch:{epoch}")

            total = 0
            correct = 0
            total_loss = 0
            # false_id = set()
            with torch.no_grad():
                total_batches = len(valid_dataloader)
                true_P = 0
                true_N = 0
                false_P = 0
                false_N = 0
                for signals, labels, reads_id in valid_dataloader:
                    signals = signals.reshape(signals.shape[0], 1, signals.shape[1])
                    signals = signals.cuda()
                    labels = labels.cuda()
                    outputs = model(signals)
                    _, predicted = torch.max(outputs.data, 1)
                    loss = criterion(outputs, labels)
                    total_loss += loss.item() * signals.size(0)
                    total += labels.size(0)
                    correct += (predicted == labels).sum().item()
                    true_P += ((predicted == 1) & (labels == 1)).sum().item()
                    true_N += ((predicted == 0) & (labels == 0)).sum().item()
                    false_P += ((predicted == 1) & (labels == 0)).sum().item()
                    false_N += ((predicted == 0) & (labels == 1)).sum().item()                 
                    # bool_tensor = labels==predicted
                    # index = torch.nonzero(bool_tensor)
                    # for i in index:
                    #     false_id.add(reads_id[i])
                    # print(f'Testing batch {batch_idx + 1}/{total_batches}; Accuracy: {correct/total}; epoch:{epoch}')

            valid_accuracy = correct / total
            average_loss = total_loss / len(valid_data)
            precision = true_P / (true_P + false_P)
            recall = true_P / (true_P + false_N)
            F1 = (2 * precision * recall) / (precision + recall)  
            print(f"Validation Accuracy: {valid_accuracy:.3f}; precision:{precision:.3f}; recall:{recall:.3f}; F1:{F1:.3f}; epoch:{epoch}")

            if valid_accuracy > max_acc and args.save:
                max_acc = valid_accuracy
                # with open('false_reads.txt', 'w') as f:
                #     for item in false_id:
                #         f.write(str(item) + '\n')
                torch.save(model.state_dict(), args.save)
        

        model.eval()
        total = 0
        correct = 0
        total_loss = 0
        with torch.no_grad():
            total_batches = len(test_dataloader)
            true_P = 0
            true_N = 0
            false_P = 0
            false_N = 0
            for signals, labels, reads_id in test_dataloader:
                signals = signals.reshape(signals.shape[0], 1, signals.shape[1])
                signals = signals.cuda()
                labels = labels.cuda()
                outputs = model(signals)
                _, predicted = torch.max(outputs.data, 1)
                loss = criterion(outputs, labels)
                total_loss += loss.item() * signals.size(0)
                total += labels.size(0)
                correct += (predicted == labels).sum().item()
                true_P += ((predicted == 1) & (labels == 1)).sum().item()
                true_N += ((predicted == 0) & (labels == 0)).sum().item()
                false_P += ((predicted == 1) & (labels == 0)).sum().item()
                false_N += ((predicted == 0) & (labels == 1)).sum().item()                 

        test_accuracy = correct / total
        average_loss = total_loss / len(test_data)
        precision = true_P / (true_P + false_P)
        recall = true_P / (true_P + false_N)
        F1 = (2 * precision * recall) / (precision + recall)  
        print(f"Test Accuracy: {test_accuracy:.3f}; precision:{precision:.3f}; recall:{recall:.3f}; F1:{F1:.3f}")


    if args.external:
        with torch.no_grad():  
            FP_reads_id = []         
            for signals, labels, reads_id in external_dataloader:
                signals = signals.reshape(signals.shape[0], 1, signals.shape[1])
                signals = signals.cuda()
                outputs = model(signals)
                _, predicted = torch.max(outputs.data, 1)
                # total += reads_id.size(0) 
                outputs = outputs.cpu()
                # print(type(reads_id))
                FP_reads_id.extend([reads_id[i] for i in range(len(reads_id)) if predicted[i]])
                # pred += [i[1] for i in outputs.data]
        FP_reads_id = set(FP_reads_id)
        if args.output:
            with open(args.output,'w') as f:
                f.write("\n".join(FP_reads_id))